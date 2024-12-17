Bragg Peaks detection using a pre-trained Faster RCNN deep neural network 
================

Inorder to use the pre-trained Faster RCNN model inside mantid using an IDAaaS instance, below steps are required.

* Launch an IDAaaS instance with GPUs from WISH > Wish Single Crystal GPU Advanced
* Launch Mantid workbench nightly from Applications->Software->Mantid->Mantid Workbench Nightly 
* Download `scriptrepository\diffraction\WISH` directory from mantid's script repository as instructed here https://docs.mantidproject.org/nightly/workbench/scriptrepository.html
* Check whether `<local path>\diffraction\WISH` path is listed under `Python Script Directories` tab from `File->Manage User Directories` of Mantid workbench.
* Below is an example code snippet to test the code. It will create a peaks workspace with the inferred peaks from the cnn. The valid values for the clustering are QLab, HDBSCAN, KMeans.
```python
from cnn.BraggDetectCNN import BraggDetectCNN
model_weights = r'/mnt/ceph/auxiliary/wish/BraggDetect_FasterRCNN_Resnet50_Weights_v1.pt'
cnn_peaks_detector = BraggDetectCNN(model_weights_path=model_weights, batch_size=64)
clustering_params = {"name":"QLab", "q_tol": 0.05}
cnn_peaks_detector.find_bragg_peaks(workspace='WISH00042730', output_ws_name="CNN_Peaks", conf_threshold=0.0, **clustering_params)
```
* If the above import is not working, check whether the `<local path>\diffraction\WISH` path is listed under `Python Script Directories` tab from `File->Manage User Directories`.
