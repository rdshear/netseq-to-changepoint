cd /Users/robertshear/Projects/netseq-to-changepoint/docker_netcpa
imagename=rdshear/netcpa
tagname=${imagename}:$(cat ../VERSION)
docker build . -t ${tagname}
docker tag ${tagname} ${imagename}:latest
docker push ${tagname}
docker push ${imagename}:latest
docker system prune -f
