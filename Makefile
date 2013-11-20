cablastp-compress: src/cablastp-compress
	cp src/cablastp-compress .

src/cablastp-compress:
	(cd src && make compress)

clean:
	(cd src && make clean)
	rm -f cablastp-*

push:
	git push origin master
	git push github master

