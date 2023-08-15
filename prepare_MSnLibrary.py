import os
from matchms.importing import load_from_mgf
#from matchms.filtering import default_filters, normalize_intensities
#from matchms import calculate_scores
#from matchms.similarity import CosineGreedy


if __name__ == '__main__':


    print('loading spectra,..')

    file_mgf = "C:/PostDoc/MSn_libraries/20230811_nih_library_pos_all_lib_MSn.mgf"
    spectrums = list(load_from_mgf(file_mgf, metadata_harmonization=False))

    print("loaded")
    print(len(spectrums))

    pass



