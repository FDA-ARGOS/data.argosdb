network:
	@docker network create -d bridge covid-drdb-network 2>/dev/null || true

builder:
	@docker build . -t hivdb/covid-drdb-builder:latest

docker-envfile:
	@test -f docker-envfile || (echo "Config file 'docker-envfile' not found, use 'docker-envfile.example' as a template to create it." && false)

update-builder:
	@docker pull hivdb/covid-drdb-builder:latest > /dev/null

inspect-builder: network docker-envfile
	@docker run --rm -it \
		--volume=$(shell pwd):/covid-drdb/ \
		--volume=$(shell dirname $$(pwd))/covid-drdb-payload:/covid-drdb/payload \
		--network=covid-drdb-network \
		--volume ~/.aws:/root/.aws:ro \
		--env-file ./docker-envfile \
   		hivdb/covid-drdb-builder:latest /bin/bash

release-builder:
	@docker push hivdb/covid-drdb-builder:latest

autofill: update-builder
	@docker run --rm \
		--volume=$(shell pwd):/covid-drdb/ \
		--volume=$(shell dirname $$(pwd))/covid-drdb-payload:/covid-drdb/payload \
   		hivdb/covid-drdb-builder:latest \
		pipenv run python -m drdb.entry autofill-payload payload/

new-study: update-builder
	@docker run --rm -it \
		--volume=$(shell pwd):/covid-drdb/ \
		--volume=$(shell dirname $$(pwd))/covid-drdb-payload:/covid-drdb/payload \
   		hivdb/covid-drdb-builder:latest \
		pipenv run python -m drdb.entry new-study payload/

import-sierra-mutations: update-builder
	@bash scripts/import_sierra_mutations.sh

import-sra-info: update-builder
	@docker run --rm -it \
		--volume=$(shell pwd):/covid-drdb/ \
		--volume=$(shell dirname $$(pwd))/covid-drdb-payload:/covid-drdb/payload \
   		hivdb/covid-drdb-builder:latest \
		pipenv run python -m drdb.entry import-sra-info payload/

sync-refaa: update-builder
	@docker run --rm \
		--volume=$(shell pwd):/covid-drdb/ \
		--volume=$(shell dirname $$(pwd))/covid-drdb-payload:/covid-drdb/payload \
   		hivdb/covid-drdb-builder:latest \
		pipenv run python -m drdb.entry update-ref-amino-acid payload/

sync-varcons: update-builder
	@docker run --rm \
		--volume=$(shell pwd):/covid-drdb/ \
		--volume=$(shell dirname $$(pwd))/covid-drdb-payload:/covid-drdb/payload \
   		hivdb/covid-drdb-builder:latest \
		pipenv run python -m drdb.entry update-variant-consensus payload/

sync-glue: update-builder
	@docker run --rm \
		--volume=$(shell pwd):/covid-drdb/ \
		--volume=$(shell dirname $$(pwd))/covid-drdb-payload:/covid-drdb/payload \
   		hivdb/covid-drdb-builder:latest \
		pipenv run python -m drdb.entry update-glue-prevalence payload/

fetch-payload-origin:
	@git -C payload/ fetch origin

local-release: update-builder network docker-envfile
	@docker run --rm \
		--shm-size=2048m \
		--volume=$(shell pwd):/covid-drdb/ \
		--volume=$(shell dirname $$(pwd))/covid-drdb-payload:/covid-drdb/payload \
		--network=covid-drdb-network \
		--volume ~/.aws:/root/.aws:ro \
		--env-file ./docker-envfile \
   		hivdb/covid-drdb-builder:latest \
		scripts/export-sqlite.sh local

release: update-builder network docker-envfile fetch-payload-origin
	@set -e ; \
	version=$$(docker run --rm -it --volume=$$(pwd):/covid-drdb/ --volume=$$(dirname $$(pwd))/covid-drdb-payload:/covid-drdb/payload --network=covid-drdb-network --volume ~/.aws:/root/.aws:ro --env-file ./docker-envfile hivdb/covid-drdb-builder:latest scripts/create-tag.sh) ; \
	git -C payload/ tag $${version} ; \
	git -C payload/ push origin $${version} ; \
	docker run \
		--rm -it \
		--volume=$$(pwd):/covid-drdb/ \
		--volume=$$(dirname $$(pwd))/covid-drdb-payload:/covid-drdb/payload \
		--env-file ./docker-envfile \
		hivdb/covid-drdb-builder:latest \
		scripts/watch-release.sh $${version}

pre-release: update-builder network docker-envfile fetch-payload-origin
	@set -e ; \
	version=$$(docker run --rm -it --volume=$$(pwd):/covid-drdb/ --volume=$$(dirname $$(pwd))/covid-drdb-payload:/covid-drdb/payload --network=covid-drdb-network --volume ~/.aws:/root/.aws:ro --env-file ./docker-envfile hivdb/covid-drdb-builder:latest scripts/create-tag.sh --prerelease) ; \
	git -C payload/ tag $${version} ; \
	git -C payload/ push origin $${version} ; \
	docker run \
		--rm -it \
		--volume=$$(pwd):/covid-drdb/ \
		--volume=$$(dirname $$(pwd))/covid-drdb-payload:/covid-drdb/payload \
		--env-file ./docker-envfile \
		hivdb/covid-drdb-builder:latest \
		scripts/watch-release.sh $${version}

devdb: update-builder network
	@docker run \
		--rm \
		--volume=$(shell pwd):/covid-drdb/ \
		--volume=$(shell dirname $$(pwd))/covid-drdb-payload:/covid-drdb/payload \
		hivdb/covid-drdb-builder:latest \
		scripts/export-sqls.sh
	$(eval volumes = $(shell docker inspect -f '{{ range .Mounts }}{{ .Name }}{{ end }}' covid-drdb-devdb))
	@mkdir -p local/sqls
	@docker rm -f covid-drdb-devdb 2>/dev/null || true
	@docker volume rm $(volumes) 2>/dev/null || true
	@docker run \
		-d --name=covid-drdb-devdb \
		-e POSTGRES_HOST_AUTH_METHOD=trust \
		-p 127.0.0.1:6543:5432 \
		--network=covid-drdb-network \
		--volume=$(shell pwd)/postgresql.conf:/etc/postgresql/postgresql.conf \
		--volume=$(shell pwd)/local/sqls:/docker-entrypoint-initdb.d \
		postgres:13.1 \
		-c 'config_file=/etc/postgresql/postgresql.conf'

log-devdb:
	@docker logs -f covid-drdb-devdb

psql-devdb:
	@docker exec -it covid-drdb-devdb psql -U postgres

psql-devdb-no-docker:
	@psql -U postgres -h localhost -p 6543

.PHONY: autofill network devdb *-devdb builder *-builder *-sqlite release pre-release debug-* sync-* update-builder new-study import-* fetch-*
