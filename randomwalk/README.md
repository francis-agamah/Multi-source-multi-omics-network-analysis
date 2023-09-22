# Random walk workflow

## Docker

### Build from Dockerfile

```
git clone https://gitlab.cmbi.umcn.nl/x-omics-twoc/randomwalk.git
cd randomwalk/
docker build -t random_walk_workflow:v0.1 . 
```


### Run docker container

Select mild, moderate or severe in run command:
``` 
docker run -v ./randomwalk/random_walk_configs:/random_walk/Data randomwalk_workflow:v0.1 'mild'|'moderate'|'severe'
```

The FULL path on your local machine to the `random_walk_configs` directory should be defined.