## This is the Park censoring directory

current: target
-include target.mk

######################################################################

Sources += censor.tex
censor.pdf: censor.tex

######################################################################

Sources += simhack.R
simhack.Rout: simhack.R

######################################################################

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff
Makefile: makestuff/Makefile
## makestuff: makestuff/Makefile
makestuff/Makefile:
	(ls ../makestuff/Makefile && /bin/ln -s ../makestuff) || git clone $(msrepo)/makestuff
	ls $@

-include makestuff/os.mk

-include makestuff/texdeps.mk
-include makestuff/wrapR.mk

-include makestuff/git.mk
-include makestuff/visual.mk
-include makestuff/projdir.mk
