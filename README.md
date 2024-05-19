# MichelsAna

Repository with Michels' electrons analysis. Aiming to make a common, configurable analyser for all geometries and purposes.

## Introduction

Run: `lar -c run_MichelsAna.fcl -s /pnfs/dune/tape_backed/dunepro/protodune-sp/full-reconstructed/2021/mc/out1/PDSPProd4a/18/80/06/50/PDSPProd4a_protoDUNE_sp_reco_stage1_p1GeV_35ms_sce_datadriven_18800650_40_20210414T012332Z.root -n1`

### Auxiliar modules or tools
I added auxiliar files to define functions and different methods so that the main module can use them in a cleaner way. Importat thing to bear in mind to do this:
- Change the CMAkeLists.txt --> `add_subdirectory(fcl)` and `SERVICE_LIBRARIES: duneana_MichelsAna; ${ART_ROOT_IO_ROOTINPUT_SOURCE}` are some examples of what you beed to add.
- `Aux.hh`:
```bash
namespace whatever
{
    class Aux
    {
        public:
            Aux(fhicl::ParameterSet const &p); // Initialize your class
            bool your_function(int i); // Initialize your functions
    }
}
```
- `Aux.cc`: must `#include "Aux.h"` and define your functions here.
```bash
namespace whatever

{
    Aux::Aux(fhicl::ParameterSet const &p){};
    bool Aux::your_function(int i)
    {
      return false;
    }
}
```
Finally in your main module you have to initialize the Aux class and use it!
```bash
YourAna::YourAna(fhicl::ParameterSet const & p) : EDAnalyzer{p}, aux_utils(new whatever::Aux(p)) //, 

...

bool result = aux_utils->your_function(5);
std::cout << "Result of laura: " << result << std::endl;
```

### Looking inside the module

* **Gen** tree:
* **Cry** tree:
* **AllReco** tree:
    0. Extract info on the PFP, tracks,showers, clusters, hits and optical flashes
    1. Loop over all the tracks
* **CutReco** tree:


## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## Authors

- [**Laura Pérez-Molina**](https://github.com/LauPM)

## License

[MIT](https://choosealicense.com/licenses/mit/)

## Acknowledgments

- [**DUNE Collaboration**](https://github.com/DUNE)
- [**Sergio Mathey Corchado**](https://github.com/mantheys)
- [**Rodrigo Álvarez Garrote**](https://github.com/rodralva)
- [**Aleena Rafique**](https://orcid.org/0000-0001-8057-4087)