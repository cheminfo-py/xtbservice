# xtb-service

![Docker Image Build and Test CI](https://github.com/cheminfo-py/xtb-service/workflows/Docker%20Image%20Build%20CI/badge.svg)

This is a webservice built using [FastAPI](https://github.com/tiangolo/fastapi).

Allows to:
- Simulate IR spectra given a SMILES string
- Get conformers

## Usage

To be usable on Heroku or Dokku, which use `PORT` environmental variables, you need to either create this environmental variable or put it into an `.env` file. For local development, the `run_docker.sh` script uses `8091`.

```bash
./build_docker # builds the docker image
./run_docker # starts the service
```

For production, you may want to use Docker compose

```bash
docker-compose up
```

### Customization

You have the option to customize the behavior of the app using environment variables:

- `IMAGINARY_FREQ_THRESHOLD`: sets the maximum energy in cm-1 for imaginary frequency (if this threshold is exceeded, the output will contain a warning)
- `MAX_ATOMS`: if the input contains more than this number of atoms, an error is thrown
- `TIMEOUT`: If the request takes longer than this time (in seconds) a `TimeOut` error is raised
- `CACHEDIR`: Sets the directory for the diskcache. It will be mounted by the docker container.

## Acknowledgments

This webservices heavily relies on [xtb](https://github.com/grimme-lab/xtb#citations) and ase. If you find this service useful and mention it in a publication, please also cite ase and xtb:
- [Ask Hjorth Larsen et al  _J. Phys.: Condens. Matter_  **2017**, 29, 273002](https://iopscience.iop.org/article/10.1088/1361-648X/aa680e/meta)
- C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme
  *WIREs Comput. Mol. Sci.*, **2020**, 11, e01493.
  DOI: [10.1002/wcms.1493](https://doi.org/10.1002/wcms.1493)
- S. Grimme, C. Bannwarth, P. Shushkov, *J. Chem. Theory Comput.*, **2017**, 13, 1989-2009.
  DOI: [10.1021/acs.jctc.7b00118](https://dx.doi.org/10.1021/acs.jctc.7b00118)
- C. Bannwarth, S. Ehlert and S. Grimme., *J. Chem. Theory Comput.*, **2019**, 15, 1652-1671.
  DOI: [10.1021/acs.jctc.8b01176](https://dx.doi.org/10.1021/acs.jctc.8b01176)
- P. Pracht, E. Caldeweyher, S. Ehlert, S. Grimme, *ChemRxiv*, **2019**, preprint.
  DOI: [10.26434/chemrxiv.8326202.v1](https://dx.doi.org/10.26434/chemrxiv.8326202.v1)
- S. Spicher and S. Grimme, *Angew. Chem. Int. Ed.*, **2020**, 59, 15665–15673
  DOI: [10.1002/anie.202004239](https://doi.org/10.1002/anie.202004239)

It also uses RDKit and the conformer search proposed by Riniker/Landrum:
-  RDKit: Open-source cheminformatics; http://www.rdkit.org
-  S. Riniker and G. A. Landrum _J. Chem. Inf. Model._ **2015**, 55, 12, 2562–257 DOI: [10.1021/acs.jcim.5b00654](https://doi.org/10.1021/acs.jcim.5b00654)
## Docs

You find docs on `http://127.0.0.1:$PORT/v1/docs.`
