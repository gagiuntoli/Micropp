Readme
======

This is a temporal readme for the tests. The tests are build automatically with cmake. 

Cmake runs the tests with ctest that do not produce outputs but reports the test
execution time and checks that the return value for the test process has
been 0. That's why inside the tests any exit failure, false assert, or abort
makes ctest to report the test as failed.

To run all the tests just:

```shell
make test
```

You can also use:

```shell
ctest
```

To run only one test just try:

```shell
ctest -R testname
```

You can see all the testnames that has been added and are available:

```shell
ctest -N
```

If you want to see the test outputs use:

```shell
ctest -R testname -V
```

There are other options that you can check in the man page for ctest.

How to add a new test:
---------------------

Simply write your test code with a **.cpp** or **.f90** extension. And add that
filename in the list you can see in CMakelists.txt in this directory.

To make your test really work as a test you should **test** different conditions
inside the code with *assert*, conditions *if + return a value !=0*, use *abort*. *assert* use is preferred because it produces a more verbose message, but remember that assert only produces code in debug mode.


