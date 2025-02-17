Bragg Peaks detection using a pre-trained Faster RCNN deep neural network 
================

Inorder to run the pre-trained Faster RCNN model via mantid inside an IDAaaS instance, below steps are required.

* Launch an IDAaaS instance with GPUs selected from WISH > Wish Single Crystal GPU Advanced
* From IDAaaS, launch Mantid workbench nightly from Applications->Software->Mantid->Mantid Workbench Nightly 
* Download `scriptrepository\diffraction\WISH` directory from mantid's script repository as instructed here https://docs.mantidproject.org/nightly/workbench/scriptrepository.html
* Check whether `<local path>\diffraction\WISH` path is listed under `Python Script Directories` tab from `File->Manage User Directories` of Mantid workbench.
* Below is an example code snippet to use the pretrained model for Bragg peak detection. It will create a peaks workspace with the inferred peaks from the model. The valid values for the `clustering` argument are `QLab` or `HDBSCAN`. For `QLab` method the default value of `q_tol=0.05` will be used for `BaseSX.remove_duplicate_peaks_by_qlab` method. 
```python
from cnn.BraggDetectCNN import BraggDetectCNN
model_weights = r'/mnt/ceph/auxiliary/wish/BraggDetect_FasterRCNN_Resnet50_Weights_v1.pt'
cnn_peaks_detector = BraggDetectCNN(model_weights_path=model_weights, batch_size=64)
cnn_peaks_detector.find_bragg_peaks(workspace='WISH00042730', output_ws_name="CNN_Peaks", conf_threshold=0.0, clustering="QLab")
```
* If the above import is not working, check whether the `<local path>\diffraction\WISH` path is listed under `Python Script Directories` tab from `File->Manage User Directories`.
* Depending on the selected `clustering` method in the above, the user can provide custom parameters using `kwargs` as shown below.
```
kwargs={"q_tol": 0.1} 
cnn_peaks_detector.find_bragg_peaks(workspace='WISH00042730', output_ws_name="CNN_Peaks", conf_threshold=0.0, clustering="QLab", **kwargs)

or 

kwargs={"cluster_selection_method": "leaf", "algorithm": "brute", "keep_ignored_labels": False} 
cnn_peaks_detector.find_bragg_peaks(workspace='WISH00042730', output_ws_name="CNN_Peaks", conf_threshold=0.0, clustering="HDBSCAN", **kwargs)
```
* The documentation for using HDBSCAN can be found here: https://scikit-learn.org/1.5/modules/generated/sklearn.cluster.HDBSCAN.html
* The documentation for using `BaseSX.remove_duplicate_peaks_by_qlab` can be found here: https://docs.mantidproject.org/nightly/techniques/ISIS_SingleCrystalDiffraction_Workflow.html
