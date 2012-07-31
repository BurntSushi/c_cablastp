cablastp-compress: src/cablastp-compress
	cp src/cablastp-compress .

src/cablastp-compress:
	(cd src && make compress)

clean:
	(cd src && make clean)

push:
	git push origin master

