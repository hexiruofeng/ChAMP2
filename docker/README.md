# ChAMP2 container

Build one image per R release; the same source tree is used for all of them.

```bash
for v in 4.2.3 4.3.3 4.4.3 4.5.3 4.6.1; do
  docker build -f docker/Dockerfile --build-arg R_VERSION=$v -t champ2:$v .
done
```

Run:

```bash
docker run --rm -it -v "$PWD:/data" champ2:4.6.1
```

Inside the container:

```r
library(ChAMP)
champ_check_env()
res <- champ.process("/data/IDAT", arraytype = "EPICv2")
```

The image installs all required packages plus the optional ones that are
available for that R/Bioconductor release; missing optional packages are
reported by `champ_check_env()` and do not prevent the core pipeline from
running.
