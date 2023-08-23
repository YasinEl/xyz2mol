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

    file_mgf = "/home/yasin/yasin/databases/20230811_nih_library_pos_all_lib_MSn.mgf"
    spectra = list(load_from_mgf(file_mgf, metadata_harmonization=False))
    print(len(spectra))

    spectra_ms3 = [spectrum for spectrum in spectra if spectrum.get("ms_level") == "3" and float(spectrum.get("precursor_purity")) >= 0.9]
    print(len(spectra_ms3))


    # matching_pairs = []

    # # Wrap the outer loop with tqdm for a progress bar
    # for i, spectrum1 in tqdm(enumerate(spectra), total=len(spectra)):
    #     for j, spectrum2 in enumerate(spectra):
    #         if i < j:
    #             precursor_mass1 = float(spectrum1.get("msn_precursor_mzs").strip('[]').split(',')[1].strip())
    #             precursor_mass2 = float(spectrum2.get("msn_precursor_mzs").strip('[]').split(',')[1].strip())

    #             collision_energy1 = float(spectrum1.get("msn_collision_energies").strip('[]').split(',')[1].strip())
    #             collision_energy2 = float(spectrum2.get("msn_collision_energies").strip('[]').split(',')[1].strip())

    #             if (within_tolerance(precursor_mass1, precursor_mass2, 20)) and (
    #                     collision_energy1 == collision_energy2):
    #                 matching_pairs.append((i, j))

    # print(f"Found {len(matching_pairs)} matching pairs")






    from concurrent.futures import ThreadPoolExecutor
    from tqdm import tqdm

    def process_spectrum(i):
        matching_pairs_local = []
        spectrum1 = spectra_ms3[i]
        for j, spectrum2 in enumerate(spectra_ms3):
            if i < j:
                precursor_mass1 = float(spectrum1.get("msn_precursor_mzs").strip('[]').split(',')[1].strip())
                precursor_mass2 = float(spectrum2.get("msn_precursor_mzs").strip('[]').split(',')[1].strip())
                
                if within_tolerance(precursor_mass1, precursor_mass2, 20):
                    collision_energy1 = float(spectrum1.get("msn_collision_energies").strip('[]').split(',')[1].strip())
                    collision_energy2 = float(spectrum2.get("msn_collision_energies").strip('[]').split(',')[1].strip())

                    if (collision_energy1 == collision_energy2):
                        matching_pairs_local.append((i, j))
        
        # Update the progress bar
        pbar.update(1)
        
        return matching_pairs_local

    matching_pairs = []

    # Define the number of worker threads
    num_workers = 60

    # Create a tqdm progress bar
    with tqdm(total=len(spectra_ms3)) as pbar:
        # Use ThreadPoolExecutor to run the tasks in parallel
        with ThreadPoolExecutor() as executor:
            for result in executor.map(process_spectrum, range(len(spectra_ms3))):
                matching_pairs.extend(result)

    print(f"Found {len(matching_pairs)} matching pairs")

    df = pd.DataFrame(matching_pairs, columns=['spectrum1_index', 'spectrum2_index'])
    df.to_csv('/home/yasin/yasin/MSn_pairs/MS3_pairs.csv', index=False)



