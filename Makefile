## This is censor

current: target
-include target.mk

# -include makestuff/perl.def

######################################################################

# Content

Sources += $(wildcard *.tex)

## censor.pdf: censor.tex

######################################################################

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff
Makefile: makestuff/Makefile
makestuff/Makefile:
	git clone $(msrepo)/makestuff
	ls $@

-include makestuff/os.mk

-include makestuff/texdeps.mk
## -include makestuff/wrapR.mk

-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk