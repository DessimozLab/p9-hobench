#!/usr/bin/env python3
"""Assemble literature/C_genomes_species.md: hand-written text (below) + generated tables from literature/raw/C/*.md"""
R="/scratch/nromashc/agentic-ai/workspaces/mywork/"
D=R+"literature/raw/C/"
panel=open(D+"panel_tables.md").read()
qfo=open(D+"qfo_table.md").read()
doi=open(D+"doi_table.md").read()
txt=open(R+"scripts/C_text_template.md").read()
txt=txt.replace("{{PANEL_TABLES}}",panel).replace("{{QFO_TABLE}}",qfo).replace("{{DOI_TABLE}}",doi)
open(R+"literature/C_genomes_species.md","w").write(txt)
print(len(txt.splitlines()),"lines")
