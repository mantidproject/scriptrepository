Bragg Peaks detection using a Faster RCNN model
================

Inorder to use the pretrained Faster RCNN model inside mantid below steps are required.

* Install mantid from conda `mamba create -n mantid_cnn -c mantid mantidworkbench`
* Activate the conda environment with `mamba activate mantid_cnn` 
* Launch workbench from `workbench` command
* Download the script repository's `scriptrepository\diffraction\WISH` directory as instructed here https://docs.mantidproject.org/nightly/workbench/scriptrepository.html 
* Close the workbench
* From command line, change the directory to the place where the scripts were downloaded ex: LocalScriptRepo\diffraction\WISH 
* Within the same conda enviroment, install pytorch dependancies by running `pip install -r requirements.txt`
* Install NVIDIA CUDA Deep Neural Network library (cuDNN) by running `conda install -c anaconda cudnn`
* Re-launch workbench from `workbench` command
* Below is an example code snippet to test the code. It will create a peaks workspace with the inferred peaks from the cnn and will do a peak filtering using the q_tol provided using `BaseSX.remove_duplicate_peaks_by_qlab`.
```python
from cnn.WISHMLBraggPeaksDetector import WISHMLBraggPeaksDetector
model_weights = r'path/to/pretrained/fasterrcnn_resnet50_model_weights.pt'
cnn_peaks_detector = WISHMLBraggPeaksDetector(model_weights_path=model_weights, batch_size=64)
cnn_peaks_detector.find_bragg_peaks(run_name='WISH00042730', q_tol=0.05)
```
