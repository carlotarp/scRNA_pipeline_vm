# Agent instructions

This folder is a config-driven scRNA-seq pipeline (breast tumours, 10x GEM-X Flex) meant to be run by a
human or by an LLM agent. Before running, annotating or changing anything, read
`pipeline/skill/scrna-pipeline/SKILL.md` and follow it: it explains the runner, the checkpoints where a
decision is needed (resolutions, hand-written annotation dictionaries, noise vs real populations) and how to
record each decision. Claude Code loads it automatically as the `scrna-pipeline` skill when
`.claude/skills/scrna-pipeline` points to it.
