import os
import numpy as np
import torch as tc
import torchvision
from torchvision import transforms
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from tqdm.notebook import tqdm, trange
import pickle

import albumentations as A
from albumentations.pytorch.transforms import ToTensorV2

# from torchvision import transforms as torchtrans  
from torchvision.models.detection.faster_rcnn import FastRCNNPredictor

# these are the helper libraries imported.
from engine import train_one_epoch, evaluate
import utils


def find_det_coordinates(spec):
    spec = spec-5
    x_before = spec//128
    y_before = spec%128
    y_new = 127-y_before
    x_new = 152*(x_before//152) + 151 - (x_before % 152)
    return int(y_new), int(x_new) # swapped to match x,y 


class WishDataSet(tc.utils.data.Dataset):
    def __init__(self, annotations_dir, transforms=None):
        self.annotations_root = annotations_dir
        self.annotation_f_names = list(sorted(os.listdir(annotations_dir)))
        self.transforms = transforms
        self.bbox_span_x_from_center = 2
        self.bbox_span_y_from_center = 10
        
    
    def __getitem__(self, idx):
        annot_f_path = os.path.join(self.annotations_root, self.annotation_f_names[idx])
        with open(annot_f_path, 'rb') as handle:
            annot_dict = pickle.load(handle)
            bin_data = np.load(annot_dict['nparr_path'])['arr_0']
            img_res = np.tile(bin_data[:,:,None], 3)
            # bin_index = annot_dict['bin_index']
            peaks = annot_dict['peaks']
            boxes = []
            labels = []
            for p in peaks:
                x,y = find_det_coordinates(p)
                # print(f"peak{p} x{x} y{y}")
                
                xmin = x - self.bbox_span_x_from_center
                ymin = y - self.bbox_span_y_from_center
                xmax = x + self.bbox_span_x_from_center
                ymax = y + self.bbox_span_y_from_center
                
                if xmin < 0:
                    xmin = 0
                
                if ymin < 0:
                    ymin = 0
                    
                if xmax > bin_data.shape[1]:
                    xmax = bin_data.shape[1]
                    
                if ymax > bin_data.shape[0]:
                    ymax = bin_data.shape[0]
                
                boxes.append([xmin, ymin, xmax, ymax])
                labels.append(1)
            
            # convert boxes into a torch.Tensor
            boxes = tc.as_tensor(boxes, dtype=tc.float32)

            # getting the areas of the boxes
            area = (boxes[:, 3] - boxes[:, 1]) * (boxes[:, 2] - boxes[:, 0])

            # suppose all instances are not crowd
            iscrowd = tc.zeros((boxes.shape[0],), dtype=tc.int64)
            labels = tc.as_tensor(labels, dtype=tc.int64)
            
            target = {}
            target["boxes"] = boxes
            target["labels"] = labels
            target["area"] = area
            target["iscrowd"] = iscrowd
            
            # image_id
            image_id = tc.tensor([idx])
            target["image_id"] = image_id
            
            if self.transforms:
                sample = self.transforms(image = img_res, bboxes = target['boxes'], labels = labels)
                img_res = sample['image']
                target['boxes'] = tc.Tensor(sample['bboxes'])
            
            # if self.transforms is not None:
            #     img, target = self.transforms(img, target)
            
            return img_res, target

    def __len__(self):
        return len(self.annotation_f_names)


stats = (np.array([1.26653515, 1.26653515, 1.26653515]), np.array([4.59509826, 4.59509826, 4.59509826]))

def get_transform(train, stats):
    if train:
        return A.Compose(
            [
                A.VerticalFlip(0.5),
                A.augmentations.geometric.resize.Resize(128, 128, p=1),
                A.Normalize(mean=stats[0], std=stats[1]),
                ToTensorV2(p=1.0)
            ],
            bbox_params={'format': 'pascal_voc', 'label_fields': ['labels']})
    else:
        return A.Compose(
            [
                A.augmentations.geometric.resize.Resize(128, 128, p=1),
                A.Normalize(mean=stats[0], std=stats[1]),
                ToTensorV2(p=1.0)
            ], 
            bbox_params={'format': 'pascal_voc', 'label_fields': ['labels']})


def get_object_detection_model(num_classes):

    # load a model pre-trained pre-trained on COCO
    model = torchvision.models.detection.fasterrcnn_resnet50_fpn(pretrained=True)
    
    # get number of input features for the classifier
    in_features = model.roi_heads.box_predictor.cls_score.in_features
    # replace the pre-trained head with a new one
    model.roi_heads.box_predictor = FastRCNNPredictor(in_features, num_classes) 

    return model

annotations_dir = "/home/wj1132075/Desktop/CNN_Model_Data/Annotations/"

dataset = WishDataSet(annotations_dir, transforms=get_transform(train=True, stats=stats))
dataset_test = WishDataSet(annotations_dir, transforms=get_transform(train=False, stats=stats))

tc.manual_seed(1)
indices = tc.randperm(len(dataset)).tolist()

# train test split
test_split = 0.2
tsize = int(len(dataset)*test_split)
dataset = tc.utils.data.Subset(dataset, indices[:-tsize])
dataset_test = tc.utils.data.Subset(dataset_test, indices[-tsize:])

# define training and validation data loaders
data_loader = tc.utils.data.DataLoader(dataset, batch_size=4, shuffle=True, num_workers=1,collate_fn=utils.collate_fn)
data_loader_test = tc.utils.data.DataLoader(dataset_test, batch_size=4, shuffle=False, num_workers=1, collate_fn=utils.collate_fn)


device = tc.device('cuda') if tc.cuda.is_available() else tc.device('cpu')

num_classes = 2

# get the model using our helper function
model = get_object_detection_model(num_classes)


model.to(device)

# construct an optimizer
params = [p for p in model.parameters() if p.requires_grad]
optimizer = tc.optim.SGD(params, lr=0.005, momentum=0.9, weight_decay=0.0005)

# and a learning rate scheduler which decreases the learning rate by
# 10x every 3 epochs
lr_scheduler = tc.optim.lr_scheduler.StepLR(optimizer, 
                                            step_size=3,
                                            gamma=0.1)# to train on gpu if selected.

num_epochs = 10

for epoch in trange(num_epochs):
    # training for one epoch
    print(f"starting epoch {epoch}")
    train_one_epoch(model, optimizer, data_loader, device, epoch, print_freq=10)
    print(f"epoch:{epoch} done")
    # update the learning rate
    lr_scheduler.step()
    print(f"{epoch} LR step done")
    # evaluate on the test dataset
    evaluate(model, data_loader_test, device=device)# training for 10 epochs
    print(f"{epoch} evalute done")
    break

