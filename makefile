main: main.o lambdaInit.o floodFill.o src/lambdaInit.h src/floodFill.h
	gcc -o main obj/main.o obj/lambdaInit.o obj/floodFill.o -lm

test_floodFill: floodFill.o test_floodFill.o
	gcc -o test_floodFill obj/test_floodFill.o obj/floodFill.o -lm

test_lambOseenShear0: floodFill.o test_lambOseenShear0.o lambdaInit.o
	gcc -o test_lambOseenShear0 obj/test_lambOseenShear0.o obj/floodFill.o obj/lambdaInit.o -lm

test_lambOseenShear1: floodFill.o test_lambOseenShear1.o lambdaInit.o
	gcc -o test_lambOseenShear1 obj/test_lambOseenShear1.o obj/floodFill.o obj/lambdaInit.o -lm

test_lambOseenShear2: floodFill.o test_lambOseenShear2.o lambdaInit.o
	gcc -o test_lambOseenShear2 obj/test_lambOseenShear2.o obj/floodFill.o obj/lambdaInit.o -lm

test_lambOseenShear3: floodFill.o test_lambOseenShear3.o lambdaInit.o
	gcc -o test_lambOseenShear3 obj/test_lambOseenShear3.o obj/floodFill.o obj/lambdaInit.o -lm

main.o: src/main.c 
	gcc -c src/main.c -o obj/main.o

lambdaInit.o: src/lambdaInit.c src/lambdaInit.h
	gcc -c src/lambdaInit.c -o obj/lambdaInit.o

floodFill.o: src/floodFill.c src/floodFill.h
	gcc -c src/floodFill.c -o obj/floodFill.o

test_floodFill.o: src/tests/test_floodFill.c src/floodFill.h
	gcc -c src/tests/test_floodFill.c -o obj/test_floodFill.o

test_lambOseenShear0.o: src/tests/test_lambOseenShear0.c src/floodFill.h src/lambdaInit.h
	gcc -c src/tests/test_lambOseenShear0.c -o obj/test_lambOseenShear0.o

test_lambOseenShear1.o: src/tests/test_lambOseenShear1.c src/floodFill.h src/lambdaInit.h
	gcc -c src/tests/test_lambOseenShear1.c -o obj/test_lambOseenShear1.o

test_lambOseenShear2.o: src/tests/test_lambOseenShear2.c src/floodFill.h src/lambdaInit.h
	gcc -c src/tests/test_lambOseenShear2.c -o obj/test_lambOseenShear2.o

test_lambOseenShear3.o: src/tests/test_lambOseenShear3.c src/floodFill.h src/lambdaInit.h
	gcc -c src/tests/test_lambOseenShear3.c -o obj/test_lambOseenShear3.o
