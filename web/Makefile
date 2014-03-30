pkgname := nexpokit
prefix := $(pkgname)
deploy := $(pkgname)

all:codes web

.phony: codes web

$(pkgname).tar.gz : ../.git/index	
	(cd ..; git archive --prefix $(prefix)/ master) | gzip > $(pkgname).tar.gz

codes : $(pkgname).tar.gz

index.md : ../README.md
	cp ../README.md index.md

web: index.md

deploy: all
	rsync -rtvp ./ ~/Dropbox/websites/purdue/codes/$(deploy)/ --exclude=Makefile --exclude=*~ --exclude=.DS_Store
