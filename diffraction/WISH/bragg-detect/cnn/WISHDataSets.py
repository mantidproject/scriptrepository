import torch as tc
import numpy as np
import albumentations as A
from albumentations.pytorch.transforms import ToTensorV2
from bragg_utils import make_3d_array
from mantid.simpleapi import Load, Rebunch, Workspace


class WISHWorkspaceDataSet(tc.utils.data.Dataset):
    def __init__(self, workspace):
        if isinstance(workspace, Workspace):
            if workspace.getAxis(0).getUnit().unitID() != "TOF":
                raise RuntimeError("Unit of the X-axis is expected to be TOF")
            ws = workspace
        elif isinstance(workspace, str):
            ws = Load(Filename=workspace, OutputWorkspace=workspace, EnableLogging=False)
        else:
            raise RuntimeError("Invalid workspace type - must be Workspace object or a name of a workspace to Load")
        self.rebunched_ws = Rebunch(InputWorkspace=ws, NBunch=3, OutputWorkspace="_cnn_rebunched", StoreInADS=False, EnableLogging=False)
        self.ws_3d = make_3d_array(self.rebunched_ws)
        print(f"Data set for {workspace} is created with shape{self.ws_3d.shape}")
        self.trans = A.Compose([A.pytorch.transforms.ToTensorV2(p=1.0)])

    def get_workspace(self):
        return self.rebunched_ws
        
    def get_ws_as_3d_array(self):
        return self.ws_3d
    
    def __len__(self):
        """
        Return number of tubes in the ws
        """
        return self.ws_3d.shape[0]

    def __getitem__(self, idx):
        tube_data = self.ws_3d[idx, ...]
        tube_data = np.tile(tube_data[:,:, None], 3).astype(np.float32)
        frame = self.trans(image=tube_data)
        return frame['image']