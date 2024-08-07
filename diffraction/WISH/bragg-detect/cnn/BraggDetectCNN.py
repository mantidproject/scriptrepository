import torchvision 
from torchvision.models.detection.faster_rcnn import FastRCNNPredictor
import torch as tc
import warnings
from cnn.WISHDataSets import WISHWorkspaceDataSet
import numpy as np
from bragg_utils import createPeaksWorkspaceFromIndices
from tqdm import tqdm
from Diffraction.single_crystal.base_sx import BaseSX
import time

class BraggDetectCNN:
    """
    Detects Bragg's peaks from a workspace using a pre trained deep learning model created using Faster R-CNN model with a ResNet-50-FPN backbone.
    Example usage is as shown below.

    #1) Set the path of the .pt file that contains the weights of the pre trained FasterRCNN model
    cnn_weights_path = r""

    # 2) Create a peaks workspace containing bragg peaks detected with a confidence greater than conf_threshold
    cnn_bragg_peaks_detector = BraggDetectCNN(model_weights_path=cnn_weights_path, batch_size=64, workers=0, iou_threshold=0.001)
    cnn_bragg_peaks_detector.find_bragg_peaks(workspace="WISH00042730", conf_threshold=0.0, q_tol=0.05)
    """

    def __init__(self, model_weights_path, batch_size=64, workers=0, iou_threshold=0.001):
        """
        :param model_weights_path: Path to the .pt file containing the weights of the pre trained CNN model
        :param batch_size: Batch size to be used when loading tube data for inferencing. 
        :param workers: Number of loader worker processes to do multi-process data loading. workers=0 means data loading in main process
        :param iou_threshold: IOU(Intersection Over Union) threshold to filter out overlapping bounding boxes for detected peaks
        """
        self.model_weights_path = model_weights_path
        self.device = self._select_device()
        self.model = self._load_cnn_model_from_weights(self.model_weights_path)
        self.batch_size = batch_size
        self.workers = workers
        self.iou_threshold = iou_threshold


    def find_bragg_peaks(self, workspace, output_ws_name="CNN_Peaks", conf_threshold=0.0, q_tol=0.05):
        """
        Find bragg peaks using the pre trained FasterRCNN model and create a peaks workspace
        :param workspace: Workspace name or the object of Workspace from WISH, ex: "WISH0042730"
        :param output_ws_name: Name of the peaks workspace
        :param conf_threshold: Confidence threshold to filter peaks inferred from RCNN
        :param q_tol: qlab tolerance to remove duplicate peaks
        """
        start_time = time.time()
        data_set, predicted_indices = self._do_cnn_inferencing(workspace)
        filtered_indices = predicted_indices[predicted_indices[:, -1] > conf_threshold]
        filtered_indices_rounded = np.round(filtered_indices[:, :-1]).astype(int)
        peaksws = createPeaksWorkspaceFromIndices(data_set.get_workspace(), output_ws_name, filtered_indices_rounded, data_set.get_ws_as_3d_array())
        for ipk, pk in enumerate(peaksws):
            pk.setIntensity(filtered_indices[ipk, -1])

        #Filter duplicates by qlab
        BaseSX.remove_duplicate_peaks_by_qlab(peaksws, q_tol)
        data_set.delete_rebunched_ws()
        print(f"Bragg peaks finding from FasterRCNN model is completed in {time.time()-start_time} seconds!")


    def _do_cnn_inferencing(self, workspace):
        data_set = WISHWorkspaceDataSet(workspace)
        data_loader = tc.utils.data.DataLoader(data_set, batch_size=self.batch_size, shuffle=False, num_workers=self.workers)
        self.model.eval()
        predicted_indices_with_score = []
        for batch_idx, img_batch in enumerate(tqdm(data_loader, desc="Processing batches of tubes")):
            for img_idx, img in enumerate(img_batch):
                tube_idx = batch_idx * data_loader.batch_size + img_idx
                with tc.no_grad():
                    prediction = self.model([img.to(self.device)])[0]
                    nms_prediction = self._apply_nms(prediction, self.iou_threshold)
                    for box, score in zip(nms_prediction['boxes'], nms_prediction['scores']):
                        tof = (box[0]+box[2])/2
                        tube_res = (box[1]+box[3])/2
                        predicted_indices_with_score.append([tube_idx, tube_res.cpu(), tof.cpu(), score.cpu()])
        return data_set, np.array(predicted_indices_with_score)
                    

    def _apply_nms(self, orig_prediction, iou_thresh):
        keep = torchvision.ops.nms(orig_prediction['boxes'], orig_prediction['scores'], iou_thresh)
        final_prediction = orig_prediction
        final_prediction['boxes'] = final_prediction['boxes'][keep]
        final_prediction['scores'] = final_prediction['scores'][keep]
        final_prediction['labels'] = final_prediction['labels'][keep]
        return final_prediction


    def _select_device(self):
        if tc.cuda.is_available():
            print("GPU device is found!")
            return tc.device("cuda")
        else:
            warnings.warn(
                "Warning! GPU is not available, the program will run very slow..", RuntimeWarning)
            return tc.device("cpu")


    def _load_cnn_model_from_weights(self, weights_path):
        model = self._get_fasterrcnn_resnet50_fpn(num_classes=2)
        model.load_state_dict(tc.load(weights_path, map_location=self.device))
        return model.to(self.device)


    def _get_fasterrcnn_resnet50_fpn(self, num_classes=2):
        model = torchvision.models.detection.fasterrcnn_resnet50_fpn(weights=None)
        in_features = model.roi_heads.box_predictor.cls_score.in_features
        # replace the head
        model.roi_heads.box_predictor = FastRCNNPredictor(in_features, num_classes) 
        return model
    