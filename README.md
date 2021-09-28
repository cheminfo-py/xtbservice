# xtb-service

![Docker Image Build and Test CI](https://github.com/cheminfo-py/xtb-service/workflows/Docker%20Image%20Build%20CI/badge.svg)

This is a webservice built using [FastAPI](https://github.com/tiangolo/fastapi).

Allows to:
- Simulate IR spectra given a SMILES string

## Usage

To be usable on Heroku or Dokku, which use `PORT`environmental variables, you need to either create this environmental variable or put it into an `.env` file. For local development, the `run_docker.sh` script uses `8091`.

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


## Docs

You find docs on `http://127.0.0.1:$PORT/docs.`
