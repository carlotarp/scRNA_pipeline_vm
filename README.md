# scRNA-seq pipeline (10x GEM-X Flex, breast tumours)

DecontX + QC, Harmony integration, lineage / leukocyte / stroma / tumour annotation, CopyKAT and
sub-population studies, driven by a config file and checkpoints where a human (or an agent) decides
resolutions and writes the cluster-to-label dictionaries.

## Install
Clone it as a folder named `pipeline` inside your project; inputs (`data/`) and outputs (`results/`)
live next to it, in the project folder:

```bash
cd my_project
git clone https://github.com/carlotarp/scRNA_pipeline_vm pipeline
mkdir -p .claude/skills && ln -s ../../pipeline/skill/scrna-pipeline .claude/skills/scrna-pipeline   # agent skill (Claude Code)
pipeline/run.sh --status
```

- How to run it, step by step: `AGENTS.md`
- Configuration: `config.yaml`, samples: `samples.csv`
- Agent skill: `skill/scrna-pipeline/`
