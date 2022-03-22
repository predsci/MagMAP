

import oftpy.utilities.datatypes.datatypes as psi_dtypes
import oftpy.utilities.plotting.psi_plotting as psi_plot

fits_path = "/Users/turtle/Dropbox/MyOFT/download_test/hmi_raw/2021/01/01/hmi_m_720s_20210101T235952_.fits"

# open fits file to a MagnetoLOS object
hmi_los = psi_dtypes.read_hmi720s(fits_file=fits_path, solar_north_up=True)

# plot the raw data
psi_plot.PlotDiskImage(disk_image=hmi_los, nfig=0, plot_attr="data", title="Raw LOS Disk")

# calculate Br
hmi_los.get_coordinates()
hmi_los.add_Br()

psi_plot.PlotDiskImage(disk_image=hmi_los, nfig=1, plot_attr="Br", title="Estimated B_{r} Disk")



