galsim: galsim.c
	gcc -O3 -ffast-math galsim.c graphics/graphics.c graphics/graphics.h -o galsim -lm -lX11
run:
	./galsim 2000 input_data/ellipse_N_02000.gal 200 0.00001 0.252 0

clean:
	rm galsim
