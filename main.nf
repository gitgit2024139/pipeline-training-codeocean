#!/usr/bin/env nextflow
// hash:sha256:849ba9c78470f84053447de51ff84e950e6da69da5f5a51efe752e0960347d8b

// capsule - FastQC
process capsule_fastqc_1 {
	tag 'capsule-4603095'
	container "registry.apps.codeocean.com/published/a77487b1-1307-43b8-8f5a-b18255789cf3:v6"

	cpus 8
	memory '7.5 GB'

	input:
	val path1

	output:
	path 'capsule/results/*.zip', emit: to_capsule_multiqc_2_2

	script:
	"""
	#!/usr/bin/env bash
	set -e

	export CO_CAPSULE_ID=a77487b1-1307-43b8-8f5a-b18255789cf3
	export CO_CPUS=8
	export CO_MEMORY=8053063680

	mkdir -p capsule
	mkdir -p capsule/data && ln -s \$PWD/capsule/data /data
	mkdir -p capsule/results && ln -s \$PWD/capsule/results /results
	mkdir -p capsule/scratch && ln -s \$PWD/capsule/scratch /scratch

	ln -s "/tmp/data/reads/$path1" "capsule/data/$path1" # id: a030c692-d697-4ed8-9c70-02ad9320a85d

	echo "[${task.tag}] cloning git repo..."
	if [[ "\$(printf '%s\n' "2.20.0" "\$(git version | awk '{print \$3}')" | sort -V | head -n1)" = "2.20.0" ]]; then
		git -c credential.helper= clone --filter=tree:0 --branch v6.0 "https://apps.codeocean.com/capsule-4603095.git" capsule-repo
	else
		git -c credential.helper= clone --branch v6.0 "https://apps.codeocean.com/capsule-4603095.git" capsule-repo
	fi
	mv capsule-repo/code capsule/code && ln -s \$PWD/capsule/code /code
	rm -rf capsule-repo

	echo "[${task.tag}] running capsule..."
	cd capsule/code
	chmod +x run
	./run ${params.capsule_fastqc_1_args}

	echo "[${task.tag}] completed!"
	"""
}

// capsule - STAR Alignment
process capsule_star_alignment_3 {
	tag 'capsule-0281418'
	container "registry.apps.codeocean.com/published/390486e3-9183-446d-a441-f61bb939c3e0:v3"

	cpus 4
	memory '30 GB'

	input:
	val path3

	output:
	path 'capsule/results/*', emit: to_capsule_sambamba_sort_and_index_4_4

	script:
	"""
	#!/usr/bin/env bash
	set -e

	export CO_CAPSULE_ID=390486e3-9183-446d-a441-f61bb939c3e0
	export CO_CPUS=4
	export CO_MEMORY=32212254720

	mkdir -p capsule
	mkdir -p capsule/data && ln -s \$PWD/capsule/data /data
	mkdir -p capsule/results && ln -s \$PWD/capsule/results /results
	mkdir -p capsule/scratch && ln -s \$PWD/capsule/scratch /scratch

	ln -s "/tmp/data/star_index" "capsule/data/star_index" # id: 3c8c7a4d-ae48-45f9-8901-14dc392ad19e
	ln -s "/tmp/data/reads/$path3" "capsule/data/$path3" # id: a030c692-d697-4ed8-9c70-02ad9320a85d

	echo "[${task.tag}] cloning git repo..."
	if [[ "\$(printf '%s\n' "2.20.0" "\$(git version | awk '{print \$3}')" | sort -V | head -n1)" = "2.20.0" ]]; then
		git -c credential.helper= clone --filter=tree:0 --branch v3.0 "https://apps.codeocean.com/capsule-0281418.git" capsule-repo
	else
		git -c credential.helper= clone --branch v3.0 "https://apps.codeocean.com/capsule-0281418.git" capsule-repo
	fi
	mv capsule-repo/code capsule/code && ln -s \$PWD/capsule/code /code
	rm -rf capsule-repo

	echo "[${task.tag}] running capsule..."
	cd capsule/code
	chmod +x run
	./run ${params.capsule_star_alignment_3_args}

	echo "[${task.tag}] completed!"
	"""
}

// capsule - MultiQC
process capsule_multiqc_2 {
	tag 'capsule-0426636'
	container "registry.apps.codeocean.com/published/5960f4da-36d4-415f-9619-c578035da1fc:v1"

	cpus 1
	memory '7.5 GB'

	publishDir "$RESULTS_PATH/quality_control", saveAs: { filename -> new File(filename).getName() }

	input:
	path 'capsule/data/'

	output:
	path 'capsule/results/*'

	script:
	"""
	#!/usr/bin/env bash
	set -e

	export CO_CAPSULE_ID=5960f4da-36d4-415f-9619-c578035da1fc
	export CO_CPUS=1
	export CO_MEMORY=8053063680

	mkdir -p capsule
	mkdir -p capsule/data && ln -s \$PWD/capsule/data /data
	mkdir -p capsule/results && ln -s \$PWD/capsule/results /results
	mkdir -p capsule/scratch && ln -s \$PWD/capsule/scratch /scratch

	echo "[${task.tag}] cloning git repo..."
	if [[ "\$(printf '%s\n' "2.20.0" "\$(git version | awk '{print \$3}')" | sort -V | head -n1)" = "2.20.0" ]]; then
		git -c credential.helper= clone --filter=tree:0 --branch v1.0 "https://apps.codeocean.com/capsule-0426636.git" capsule-repo
	else
		git -c credential.helper= clone --branch v1.0 "https://apps.codeocean.com/capsule-0426636.git" capsule-repo
	fi
	mv capsule-repo/code capsule/code && ln -s \$PWD/capsule/code /code
	rm -rf capsule-repo

	echo "[${task.tag}] running capsule..."
	cd capsule/code
	chmod +x run
	./run ${params.capsule_multiqc_2_args}

	echo "[${task.tag}] completed!"
	"""
}

// capsule - Sambamba Sort and Index
process capsule_sambamba_sort_and_index_4 {
	tag 'capsule-0252197'
	container "registry.apps.codeocean.com/published/b801a697-4549-4874-83f1-191f39f37ae4:v1"

	cpus 1
	memory '7.5 GB'

	input:
	path 'capsule/data/'

	output:
	path 'capsule/results/*', emit: to_capsule_featurecounts_5_5

	script:
	"""
	#!/usr/bin/env bash
	set -e

	export CO_CAPSULE_ID=b801a697-4549-4874-83f1-191f39f37ae4
	export CO_CPUS=1
	export CO_MEMORY=8053063680

	mkdir -p capsule
	mkdir -p capsule/data && ln -s \$PWD/capsule/data /data
	mkdir -p capsule/results && ln -s \$PWD/capsule/results /results
	mkdir -p capsule/scratch && ln -s \$PWD/capsule/scratch /scratch

	echo "[${task.tag}] cloning git repo..."
	if [[ "\$(printf '%s\n' "2.20.0" "\$(git version | awk '{print \$3}')" | sort -V | head -n1)" = "2.20.0" ]]; then
		git -c credential.helper= clone --filter=tree:0 --branch v1.0 "https://apps.codeocean.com/capsule-0252197.git" capsule-repo
	else
		git -c credential.helper= clone --branch v1.0 "https://apps.codeocean.com/capsule-0252197.git" capsule-repo
	fi
	mv capsule-repo/code capsule/code && ln -s \$PWD/capsule/code /code
	rm -rf capsule-repo

	echo "[${task.tag}] running capsule..."
	cd capsule/code
	chmod +x run
	./run ${params.capsule_sambamba_sort_and_index_4_args}

	echo "[${task.tag}] completed!"
	"""
}

// capsule - featureCounts
process capsule_featurecounts_5 {
	tag 'capsule-1555597'
	container "registry.apps.codeocean.com/published/aa8a9cb9-893e-4d0b-8194-bd88fee2fe62:v2"

	cpus 2
	memory '15 GB'

	publishDir "$RESULTS_PATH/feature_counts/1", saveAs: { filename -> new File(filename).getName() }

	input:
	path 'capsule/data/'

	output:
	path 'capsule/results/*'

	script:
	"""
	#!/usr/bin/env bash
	set -e

	export CO_CAPSULE_ID=aa8a9cb9-893e-4d0b-8194-bd88fee2fe62
	export CO_CPUS=2
	export CO_MEMORY=16106127360

	mkdir -p capsule
	mkdir -p capsule/data && ln -s \$PWD/capsule/data /data
	mkdir -p capsule/results && ln -s \$PWD/capsule/results /results
	mkdir -p capsule/scratch && ln -s \$PWD/capsule/scratch /scratch

	ln -s "/tmp/data/annotation" "capsule/data/annotation" # id: 16326904-0e06-4086-841f-a4b0999d88ce

	echo "[${task.tag}] cloning git repo..."
	if [[ "\$(printf '%s\n' "2.20.0" "\$(git version | awk '{print \$3}')" | sort -V | head -n1)" = "2.20.0" ]]; then
		git -c credential.helper= clone --filter=tree:0 --branch v2.0 "https://apps.codeocean.com/capsule-1555597.git" capsule-repo
	else
		git -c credential.helper= clone --branch v2.0 "https://apps.codeocean.com/capsule-1555597.git" capsule-repo
	fi
	mv capsule-repo/code capsule/code && ln -s \$PWD/capsule/code /code
	rm -rf capsule-repo

	echo "[${task.tag}] running capsule..."
	cd capsule/code
	chmod +x run
	./run ${params.capsule_featurecounts_5_args}

	echo "[${task.tag}] completed!"
	"""
}

workflow {
	// input data
	reads_to_fastqc_1 = Channel.fromPath("../data/reads/*", type: 'any', relative: true)
	reads_to_star_alignment_3 = Channel.fromPath("../data/reads/*", type: 'any', relative: true)

	// run processes
	capsule_fastqc_1(reads_to_fastqc_1)
	capsule_star_alignment_3(reads_to_star_alignment_3)
	capsule_multiqc_2(capsule_fastqc_1.out.to_capsule_multiqc_2_2.collect())
	capsule_sambamba_sort_and_index_4(capsule_star_alignment_3.out.to_capsule_sambamba_sort_and_index_4_4)
	capsule_featurecounts_5(capsule_sambamba_sort_and_index_4.out.to_capsule_featurecounts_5_5.collect())
}
