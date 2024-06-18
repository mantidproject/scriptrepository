import torch as tc
import os
import pickle
import numpy as np
import albumentations as A
from albumentations.pytorch.transforms import ToTensorV2
from wish_utils import *

def get_transform(train):
    if train:
        return A.Compose(
            [
                A.VerticalFlip(0.5),
                # A.HorizontalFlip(0.5),
                ToTensorV2(p=1.0)
            ],
            bbox_params={'format': 'pascal_voc', 'label_fields': ['labels']})
    else:
        return A.Compose(
            [
                ToTensorV2(p=1.0)
            ], 
            bbox_params={'format': 'pascal_voc', 'label_fields': ['labels']})



class WishWorkspaceDataSet(tc.utils.data.Dataset):
    def __init__(self, ws_path):
        self.ws_path = ws_path
        self.ws = np.load(self.ws_path)['arr_0']
        print(f"data set create for {self.ws_path} shape{self.ws.shape}")
        self.trans = A.Compose([ToTensorV2(p=1.0)])

    def __len__(self):
        """
        Return number of tubes in the ws
        """
        return self.ws.shape[0]

    def __getitem__(self, idx):
        tube_data = self.ws[idx, ...]
        tube_data = np.tile(tube_data[:,:, None], 3).astype(np.float32)
        frame = self.trans(image=tube_data)
        return frame['image']
    
    


##################################################################################
class WishTubeFrameDataSet(tc.utils.data.Dataset):
    def __init__(self, annotations_dir, transforms=None):
        self.annotations_root = annotations_dir
        self.annotation_f_names = list(sorted(os.listdir(annotations_dir)))
        self.transforms = transforms
        self.bbox_span_y_from_center = 5

    def __getitem__(self, idx):
        annot_f_path = os.path.join(self.annotations_root, self.annotation_f_names[idx])
        with open(annot_f_path, 'rb') as handle:
            annot_dict = pickle.load(handle)
            frame_data = np.load(annot_dict['nparr_path'])['arr_0']
            img_res = np.tile(frame_data[:,:,None], 3).astype(np.float32)
            
            peaks = annot_dict['peaks']
            boxes = []
            labels = []
            for y, s, x, e  in peaks:
                xmin = s
                ymin = y - self.bbox_span_y_from_center
                xmax = e
                ymax = y + self.bbox_span_y_from_center
                
                if xmin < 0:
                    xmin = 0
                
                if ymin < 0:
                    ymin = 0
                    
                if xmax > frame_data.shape[1]:
                    xmax = frame_data.shape[1]
                    
                if ymax > frame_data.shape[0]:
                    ymax = frame_data.shape[0]
                
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
            target["image_id"] = idx
            
            if self.transforms:
                sample = self.transforms(image = img_res, bboxes = target['boxes'], labels = labels)
                img_res = sample['image']
                target['boxes'] = tc.Tensor(sample['bboxes'])
            
            return img_res, target

    def __len__(self):
        return len(self.annotation_f_names)
        

###########################################################################################
class WishDataSet(tc.utils.data.Dataset):
    def __init__(self, annotations_dir, transforms=None):
        self.annotations_root = annotations_dir
        self.annotation_f_names = list(sorted(os.listdir(annotations_dir)))
        self.transforms = transforms
        self.bbox_span_x_from_center = 5
        self.bbox_span_y_from_center = 3
        
    
    def __getitem__(self, idx):
        annot_f_path = os.path.join(self.annotations_root, self.annotation_f_names[idx])
        with open(annot_f_path, 'rb') as handle:
            annot_dict = pickle.load(handle)
            bin_data = np.load(annot_dict['nparr_path'])['arr_0']
            img_res = np.tile(bin_data[:,:,None], 3).astype(np.float32)
            
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
            
            # image_id = tc.tensor([idx])
            target["image_id"] = idx
            
            if self.transforms:
                sample = self.transforms(image = img_res, bboxes = target['boxes'], labels = labels)
                img_res = sample['image']
                target['boxes'] = tc.Tensor(sample['bboxes'])
            
            return img_res, target

    def __len__(self):
        return len(self.annotation_f_names)