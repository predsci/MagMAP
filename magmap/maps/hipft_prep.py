
"""
These functions are used for taking a raw Br map and preparing it for assimilation
into HIPFT.  In the simplest case, this means assigning an assimilation-weight to each
pixel.  Assimilation weights vary from 0 (do not use) to 1 (new data overwrites model
state).
Programmatically, this stage of data prep consists of layer creation and management in
the MagnetoMap object.  We always expect the first layer to be Br, and the second to
be assimilation-weight.  Additional layers can be added, but must be registered in
MagnetoMap.layers in order to be passed to HIPFT.  Layers are stored in the MagnetoMap
object as individual 2D numpy.ndarray, but will be passed to HIPFT as a single 3D
half-precision array.  The order of layers is determined by MagnetoMap.layers.
"""
import sys

import pandas as pd
import numpy as np


def set_assim_wghts(br_map, assim_method="mu4_upton"):
    # wrapper for various data-side assimilation methods
    if assim_method == "mu":
        # use mu as the assimilation weight
        br_map.assim_weight = br_map.mu
        # give pixels from the backside of the sun a 0 weight
        br_map.assim_weight[br_map.mu < 0.] = 0.
        # give pixels with 'no_data_val' a 0 weight
        br_map.assim_weight[br_map.data == br_map.no_data_val] = 0
        # update layer order
        br_map.layers = br_map.layers.append(pd.DataFrame(dict(layer_num=[1, 2], layer_name=["assim_weight", "mu"],
                                                               py_field_name=['assim_weight', "mu"])))
    elif assim_method == "mu4":
        # use mu^4 as the assimilation weight
        br_map.assim_weight = np.power(br_map.mu, 4)
        # give pixels from the backside of the sun a 0 weight
        br_map.assim_weight[br_map.mu < 0.] = 0.
        # give pixels with 'no_data_val' a 0 weight
        br_map.assim_weight[br_map.data == br_map.no_data_val] = 0.
        # update layer order
        br_map.layers = br_map.layers.append(pd.DataFrame(dict(layer_num=[1, 2], layer_name=["assim_weight", "mu"],
                                                               py_field_name=['assim_weight', "mu"])))
    elif assim_method == "mu4_upton":
        # use mu^4 as the assimilation weight (mu <= 0.1 gets no weight)
        br_map.assim_weight = np.power(br_map.mu, 4)
        # give pixels mu <= 0.1 a 0 weight (and the backside of the sun)
        br_map.assim_weight[br_map.mu <= 0.1] = 0.
        # give pixels with 'no_data_val' a 0 weight
        br_map.assim_weight[br_map.data == br_map.no_data_val] = 0
        # update layer order
        br_map.layers = br_map.layers.append(pd.DataFrame(dict(layer_num=[1, 2], layer_name=["assim_weight", "mu"],
                                                               py_field_name=['assim_weight', "mu"])))
    elif assim_method == "custom_example":
        # This example shows how to add custom layers to the map object.
        # As long as these layers are added to br_map.layers, they will be written to the
        # HIPFT .h5 file

        # use some function of mu as the assimilation weight
        br_map.assim_weight = np.power(br_map.mu, 3)
        # give pixels from the backside of the sun a 0 weight
        br_map.assim_weight[br_map.mu < 0.] = 0.
        # give pixels with 'no_data_val' a 0 weight
        br_map.assim_weight[br_map.data == br_map.no_data_val] = 0
        # update layer order
        br_map.layers = br_map.layers.append(pd.DataFrame(dict(layer_num=[1, ], layer_name=["assim_weight", ],
                                                               py_field_name=['assim_weight', ])))
        # also calculate a measurement uncertainty
        # br_map.MyLayerName must match with py_field_name=['MyLayerName', ]
        br_map.uncert = .03*br_map.data
        br_map.uncert[br_map.data == br_map.no_data_val] = br_map.no_data_val
        # update layer order
        br_map.layers = br_map.layers.append(pd.DataFrame(dict(layer_num=[2, ], layer_name=["uncertainty", ],
                                                               py_field_name=['uncert', ])))

        # pass mu as an additional layer (already exists in br_map) by adding it to the layer order
        br_map.layers = br_map.layers.append(pd.DataFrame(dict(layer_num=[3, ], layer_name=["mu", ],
                                                               py_field_name=['mu', ])))
    else:
        sys.exit("Assimilation method '" + assim_method + "' is not supported.")
    # record assimilation method in map object
    br_map.assim_method = assim_method

    return br_map



