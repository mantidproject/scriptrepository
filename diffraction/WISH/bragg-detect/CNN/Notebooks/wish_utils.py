import numpy as np
import torch as tc
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import torchvision
from torchvision.models.detection.faster_rcnn import FastRCNNPredictor, GeneralizedRCNNTransform
from tqdm.notebook import tqdm

def plot_bin(bin_data, boxes, bin_id):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    cs = plt.contourf(bin_data)
    plt.colorbar(cs)
    for box in boxes:
        det_x, det_y = (box[0]+box[2])/2, (box[1]+box[3])/2
        # plt.plot(det_y, det_x, marker="x", markersize=5, markeredgecolor="red", markerfacecolor="red")
        x, y, width, height  = box[0], box[1], box[2]-box[0], box[3]-box[1]
        rect = patches.Rectangle((x, y),
                                 width, height,
                                 linewidth = 1,
                                 edgecolor = 'r',
                                 facecolor = 'none')

        # Draw the bounding box on top of the image
        ax.add_patch(rect)
    plt.title(f"TOF bin:{bin_id}")
    plt.show()


def find_det_coordinates(spec):
    spec = spec-5
    x_before = spec//128
    y_before = spec%128
    y_new = 127-y_before
    x_new = 152*(x_before//152) + 151 - (x_before % 152)
    return int(y_new), int(x_new) # swapped to match x,y 


def get_object_detection_model(num_classes, pretrained=True):
    # load a model pre-trained pre-trained on COCO
    model = torchvision.models.detection.fasterrcnn_resnet50_fpn(pretrained=pretrained)
    # get number of input features for the classifier
    in_features = model.roi_heads.box_predictor.cls_score.in_features
    # replace the pre-trained head with a new one
    model.roi_heads.box_predictor = FastRCNNPredictor(in_features, num_classes) 
    return model


def get_object_detection_model_v2(num_classes, pretrained=True):
    model = torchvision.models.detection.fasterrcnn_resnet50_fpn_v2(pretrained=pretrained)
    in_features = model.roi_heads.box_predictor.cls_score.in_features
    model.roi_heads.box_predictor = FastRCNNPredictor(in_features, num_classes) 
    return model
    

def get_model(model, num_classes, stats=None):
    in_features = model.roi_heads.box_predictor.cls_score.in_features
    model.roi_heads.box_predictor = FastRCNNPredictor(in_features, num_classes)

    if stats != None:
        model.transform = GeneralizedRCNNTransform(min_size=800, max_size=1333, mode='bilinear', image_mean=stats[0], image_std=stats[1])

    return model


def get_stats(dl):
    cnt, csum, csum_sq = np.zeros(1), np.zeros(3), np.zeros(3)  # accumulators; 3d because 3 colour channels
    
    for xb,_ in tqdm(dl):
        xb = np.array(xb)
        cnt += xb.shape[0]*xb.shape[1]*xb.shape[2]  # number of pixels per channel
        csum += np.sum(xb, axis=(0,1,2))            # this gives 3 numbers i.e. sums of pixel intensities per each colour channel
        csum_sq += np.sum(xb**2, axis=(0,1,2))      

    μ = csum / cnt
    σ = np.sqrt(csum_sq / cnt - (μ**2))
    return μ, σ


def apply_nms(orig_prediction, iou_thresh=0.3):
    # torchvision returns the indices of the bboxes to keep
    keep = torchvision.ops.nms(orig_prediction['boxes'], orig_prediction['scores'], iou_thresh)
    
    final_prediction = orig_prediction
    final_prediction['boxes'] = final_prediction['boxes'][keep]
    final_prediction['scores'] = final_prediction['scores'][keep]
    final_prediction['labels'] = final_prediction['labels'][keep]
    
    return final_prediction

# function to convert a torchtensor back to PIL image
def torch_to_pil(img):
    return torchtrans.ToPILImage()(img).convert('RGB')


def plot_bin_detections(bin_data, real_boxes, predicted_boxes, score_threshold = 0.75):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    cs = plt.contourf(bin_data)
    plt.colorbar(cs)
    for box in real_boxes:
        det_x, det_y = (box[0]+box[2])/2, (box[1]+box[3])/2
        plt.plot(det_x, det_y, marker="x", markersize=9, markeredgecolor="lime", markerfacecolor="lime")
        # x, y, width, height  = box[0], box[1], box[2]-box[0], box[3]-box[1]
        # rect = patches.Rectangle((x, y),
        #                          width, height,
        #                          linewidth = 2,
        #                          edgecolor = 'r',
        #                          facecolor = 'none')
        # # Draw the bounding box on top of the image
        # ax.add_patch(rect)
        

    plotted_predictions = 0
    for box, score in zip(predicted_boxes['boxes'], predicted_boxes['scores']):
        print(f"prediction={box} score={score}")
        x_hat, y_hat = (box[0]+box[2])/2, (box[1]+box[3])/2
        if score < score_threshold:
            print(f"----->Ignoring the prediction at x={x_hat} y={y_hat} with score={score}")
            plt.plot(x_hat, y_hat, marker="x", markersize=5, markeredgecolor="white", markerfacecolor="white")
            continue
            
        print(f"Prediction coordinate x={x_hat} y={y_hat}")
        
        # x, y, width, height  = box[0], box[1], box[2]-box[0], box[3]-box[1]
        # rect = patches.Rectangle((x, y),
        #                          width, height,
        #                          linewidth = 1,
        #                          edgecolor = 'y',
        #                          facecolor = 'none')

        # Draw the bounding box on top of the image
        # ax.add_patch(rect)
        plt.plot(x_hat, y_hat, marker="x", markersize=5, markeredgecolor="red", markerfacecolor="red")
        plotted_predictions += 1

    plt.title(f"Predicted-{plotted_predictions} vs Real-{len(real_boxes)}")
    plt.show()


def plot_all_detections(bin_data, real_boxes, predicted_boxes):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    cs = plt.contourf(bin_data)
    plt.colorbar(cs)
    for box in real_boxes:
        det_x, det_y = (box[0]+box[2])/2, (box[1]+box[3])/2
        # x, y, width, height  = box[0], box[1], box[2]-box[0], box[3]-box[1]
        # rect = patches.Rectangle((x, y),
        #                          width, height,
        #                          linewidth = 2,
        #                          edgecolor = 'r',
        #                          facecolor = 'none')
        # ax.add_patch(rect)
        plt.plot(det_x, det_y, marker="x", markersize=9, markeredgecolor="lime", markerfacecolor="lime")

    plotted_predictions = 0
    for box, score in zip(predicted_boxes['boxes'], predicted_boxes['scores']):
        print(f"prediction={box} score={score}")
        x_hat, y_hat = (box[0]+box[2])/2, (box[1]+box[3])/2        
        print(f"Prediction coordinate x={x_hat} y={y_hat}")
        # x, y, width, height  = box[0], box[1], box[2]-box[0], box[3]-box[1]
        # rect = patches.Rectangle((x, y),
        #                          width, height,
        #                          linewidth = 1,
        #                          edgecolor = 'y',
        #                          facecolor = 'none')
        # # Draw the bounding box on top of the image
        # ax.add_patch(rect)
        plt.plot(x_hat, y_hat, marker="x", markersize=5, markeredgecolor="red", markerfacecolor="red")
        plotted_predictions += 1

    plt.title(f"Predicted-{plotted_predictions} vs Real-{len(real_boxes)}")
    plt.show()


def freeze(md, fr=True):
    ch = list(md.children())
    for c in ch: freeze(c, fr)
    # not freezing the BatchNorm layers!
    if not ch and (not isinstance(md, tc.nn.modules.batchnorm.BatchNorm2d)) and (not isinstance(md, torchvision.ops.misc.FrozenBatchNorm2d)):  
        for p in md.parameters(): 
            p.requires_grad = not fr
            # print('---\n', md, p.requires_grad)
			
			
def freeze_to(md, ix=-1, fr=True):
    ch_all = list(md.children())
    for ch in ch_all[:ix]:
        freeze(ch, fr)


def calculate_confusion_matrix_per_bin(bin_data, real_bboxes, predicted_bboxes, x_tol=2, y_tol=2):
    real_peaks = []
    predicted_peaks = []
    for box in real_bboxes:
        x = tc.round((box[0]+box[2])/2)
        y = tc.round((box[1]+box[3])/2)
        real_peaks.append((x,y))

    for box in predicted_bboxes:
        x = tc.round((box[0]+box[2])/2)
        y = tc.round((box[1]+box[3])/2)
        predicted_peaks.append((x,y))

    true_positives = 0
    false_positives = 0
    false_negatives = 0
    true_negatives = 0

    for peak in predicted_peaks:
        for r_peak in real_peaks:
            x, y = peak
            r_x, r_y = r_peak
            if abs(x-r_x) <= x_tol and abs(y-r_y) <= y_tol:
                true_positives += 1
            else:
                false_positives += 1

    for r_peak in real_peaks:
        for peak in predicted_peaks:
            x, y = peak
            r_x, r_y = r_peak
            if not(abs(x-r_x) <= x_tol and abs(y-r_y) <= y_tol):
                false_negatives += 1

    N = bin_data.shape[0] * bin_data.shape[1]
    true_negatives = N - true_positives - false_positives - false_negatives
    return true_positives, false_positives, false_negatives, true_negatives


def get_false_positive_rate(fp, tn):
    print("False Positive Rate")
    return fp/(fp+tn)*100


def get_true_positive_rate(tp, fn):
    print("True Positive Rate")
    return tp/(tp+fn)*100


def get_precision(tp, fp):
    print("Precision")
    return tp/(tp+fp)*100


def get_recall(tp, fn):
    print("Recall")
    return tp/(tp+fn)*100