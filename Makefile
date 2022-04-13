test:
	nosetests --with-doctest

increment-minor:
	bump2version minor --allow-dirty

increment-major:
	bump2version minor --allow-dirty

increment-patch:
	bump2version patch --allow-dirty