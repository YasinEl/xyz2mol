import os
from matchms.importing import load_from_mgf
from tqdm import tqdm
#from matchms.filtering import default_filters, normalize_intensities
#from matchms import calculate_scores
#from matchms.similarity import CosineGreedy


def within_tolerance(mass1, mass2, ppm):
    return abs(mass1 - mass2) <= (ppm / 1e6) * 0.5 * (mass1 + mass2)

if __name__ == '__main__':


    print('loading spectra,..')

    file_mgf = "C:/PostDoc/MSn_libraries/20230811_nih_library_pos_all_lib_MSn.mgf"
    spectra = list(load_from_mgf(file_mgf, metadata_harmonization=False))
    print(len(spectra))

    spectra = [spectrum for spectrum in spectra if spectrum.get("mslevel") == "3"]
    print(len(spectra))


    matching_pairs = []

    # Wrap the outer loop with tqdm for a progress bar
    for i, spectrum1 in tqdm(enumerate(spectra), total=len(spectra)):
        for j, spectrum2 in enumerate(spectra):
            if i < j:
                precursor_mass1 = float(spectrum1.get("msn_precursor_mzs").strip('[]').split(',')[1].strip())
                precursor_mass2 = float(spectrum2.get("msn_precursor_mzs").strip('[]').split(',')[1].strip())

                collision_energy1 = float(spectrum1.get("msn_collision_energies").strip('[]').split(',')[1].strip())
                collision_energy2 = float(spectrum2.get("msn_collision_energies").strip('[]').split(',')[1].strip())

                if (within_tolerance(precursor_mass1, precursor_mass2, 20)) and (
                        collision_energy1 == collision_energy2):
                    matching_pairs.append((i, j))

    print(f"Found {len(matching_pairs)} matching pairs")
