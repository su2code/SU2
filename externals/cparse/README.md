# CParser • [![Build Status][travis-image]][travis] [![License][license-image]][license]

[travis-image]: https://travis-ci.org/cparse/cparse.png?branch=master
[travis]: http://travis-ci.org/cparse/cparse

[release-image]: http://img.shields.io/badge/release-0.2.1-blue.svg?style=flat
[releases]: https://github.com/cmusatyalab/openface/releases

[license-image]: http://img.shields.io/badge/license-MIT-blue.svg?style=flat
[license]: LICENSE.mit

This project provides a C++ library to parse a character sequence
as an expression using Dijkstra's
[Shunting-yard algorithm](http://en.wikipedia.org/wiki/Shunting-yard_algorithm),
which modifies
[Jesse Brown's original code](http://www.daniweb.com/software-development/cpp/code/427500/calculator-using-shunting-yard-algorithm).

*This project was developed by [Brandon Amos](http://bamos.github.io) and Vinícius Garcia.*

# Getting Started

If you want to use this library in your project please take a look at our [Wiki][wiki]

[wiki]: https://github.com/cparse/cparse/wiki

## Builtin Features
 + Unary operators. +, -
 + Binary operators. +, -, /, *, %, <<, >>, ^
 + Boolean operators. <, >, <=, >=, ==, !=, &&, ||
 + Functions. sin, cos, tan, abs, print
 + Support for an hierarchy of scopes with local scope, global scope etc.
 + Easy to add new operators, operations, functions and even new types
 + Easy to implement object-to-object inheritance (with the prototype concept)
 + Built-in garbage collector (does not handle cyclic references yet)


## Setup

### Download and Compile

```bash
cd 'my/project/dir'
git clone https://github.com/cparse/cparse.git
make release -C cparse
```

### Link with your project:

```bash
g++ cparse/builtin-features.o cparse/core-shunting-yard.o main.cpp -o main
```

## Customizing your Library
To customize your calculator:

 1. Copy the `builtin-features.cpp` file and `builtin-features/` directory to your project.
 2. Edit the `builtin-features/*.inc` files as you like.
 3. Then build the project:
    1. Compile the library: `make release -C cparse/`
    2. Compile your modified features: `g++ -c builtin-features.cpp -o my-features.o`
    3. Link your project: `g++ my-features.o cparse/core-shunting-yard.o main.cpp -o main`

For a more detailed guide read our [Wiki][wiki] advanced concepts' section:
 
 + [Defining New Functions](https://github.com/bamos/cpp-expression-parser/wiki/Defining-New-Functions)
 + [Defining New Operations](https://github.com/bamos/cpp-expression-parser/wiki/Defining-New-Operations)
 + [Defining New Reserved Words](https://github.com/bamos/cpp-expression-parser/wiki/Defining-Reserved-Words)
 
## Minimal examples

### As a simple calculator

```C++
#include <iostream>
#include "shunting-yard.h"

int main() {
  TokenMap vars;
  vars["pi"] = 3.14;
  std::cout << calculator::calculate("-pi+1", &vars) << std::endl;

  // Or if you want to evaluate an expression
  // several times efficiently:
  calculator c1("pi-b");
  vars["b"] = 0.14;
  std::cout << c1.eval(vars) << std::endl; // 3
  vars["b"] = 2.14;
  std::cout << c1.eval(vars) << std::endl; // 1

  return 0;
}
```

### As a sub-parser for a programming language

Here we implement an interpreter for multiple expressions, the delimiter used
will be `;` or `\n` just like Javascript or Python and the code must start and end on curly brackets.

A similar architecture can be used for interpreting other common programming language statements like `for` loops and `if` statements. If you're interested take a look on the [jSpy programming language](https://github.com/vingarcia/jspy) that uses this project as the core parsing system.

```C++
#include <iostream>
#include "shunting-yard.h"
#include "shunting-yard-exceptions.h"

struct codeBlock {
  static void interpret(const char* start, const char** end, TokenMap vars) {
    // Remove white spaces:
    while (isspace(*start)) ++start;

    if (*start != '{') {
      throw syntax_error("Expected '{'");
    } else {
      ++start;
    }

    while (*start != '}') {
      calculator::calculate(start, vars, ";\n}", &start);

      // Alternatively you could write above:
      // - calculator(start, ";\n}", &start).eval(vars);

      // Find the beginning of the next expression:
      while(isspace(*start) || *start == ';') ++start;
    }

    if (*start == '}') {
      *end = start+1;
    } else {
      throw syntax_error("Expected '}'");
    }
  }
};

int main() {
  GlobalScope vars;
  const char* code =
    "{"
    "  a = 10;"
    "  b = 20\n"
    "  c = a + b }";

  codeBlock::interpret(code, &code, vars);

  std::cout << vars["c"] << std::endl; // 30
  return 0;
}
```

Please note that a calculator can compile an expression so that it can efficiently be executed several times at a later moment.

## More examples

 + For more examples and a comprehensible guide please read our [Wiki][wiki]
 
## Contributing

- I would like to keep this library minimal so new features should be very useful to be accepted.
- If proposed change is not a common use case, I will probably not accept it.
