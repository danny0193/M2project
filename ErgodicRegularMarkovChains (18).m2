


newPackage(
    "ErgodicRegularMarkovChains",
    Version => "0.1",
    Date => "Fall 2021",
    Headline => "Ergodic and Regular Markov Chains",
    Authors => {{ Name => "Danny Lin", Email => "dlin0710@lions.piedmont.edu", HomePage => ""}}
    )

export {"isRegularMarkov", "findingw", "createW", "TransitionMatrix", "transitionMatrix"}

-* Code section *-
TransitionMatrix = new SelfInitializingType of Matrix

transitionMatrix = method()
transitionMatrix Matrix := A -> (
    if numRows A!= numColumns A then
       error "A square matrix is expected.";
    if not all(sum entries transpose A, x -> x==1) then
       error "The rows don't add up to 1.";
    TransitionMatrix A
    )

isRegularMarkov = method()

--normal test for regular
isRegularMarkov(TransitionMatrix, ZZ) := (P,n) -> (
    I := id_(target P);
    if rank kernel(P-I)!=1 then return false
    else for i from 1 to n do(
        if not any(flatten entries P^i, zero)
        then return true);
    error"Can't be determined"
    )

isRegularMarkov(Matrix, ZZ) := (P, n) -> isRegularMarkov(TransitionMatrix P, n)

--test for regular if n isn't given
isRegularMarkov(Matrix) := (P) -> (
    output := isRegularMarkov(TransitionMatrix P)
    )

isRegularMarkov(TransitionMatrix) := (P) -> isRegularMarkov(P, 10)

findingw = method()

--outputs the nx1 matrix that each row of the transition matrix approaches
findingw(TransitionMatrix) := (P) -> (
    I := id_(target P);
    NP := transpose(P-I);
    Null := transpose gens kernel NP;
    S := sum flatten entries(Null);
    NNull := (S^-1)*(Null)
    )

findingw(Matrix) := (P) -> findingw(TransitionMatrix P)

createW = method()

--outputs the nxn matrix with identical rows
createW(TransitionMatrix) := (P) -> (
    w := findingw P;
    n := numRows P;
    W := matrix toList(n:flatten entries w)
    )

createW(Matrix) := (P) -> createW(TransitionMatrix P)


-* Documentation section *-
beginDocumentation()

doc ///
  Key
    isRegularMarkov
    (isRegularMarkov, TransitionMatrix, ZZ)
    (isRegularMarkov, Matrix, ZZ)
    (isRegularMarkov, TransitionMatrix)
    (isRegularMarkov, Matrix)
  Usage
    isRegularMarkov(P,n)
  Description
    Text
      This function requires the input of P (the transition matrix) and an optional n
      (a positive integer that is the upper bound we want to test this property
      for) and returns True for Regular Markov Chains, and False for Non Regular
      Markov Chains.
    Example
      P := matrix{{.5, .25, .25}, {.5, 0, .5}, {.25, .25, .5}}
      isRegularMarkov(P, 15)
        
///

doc ///
  Key
    findingw
    (findingw, TransitionMatrix)
    (findingw, Matrix)
  Usage
    findingw(P)
  Description
    Text
      This function requires the input of P (the transition matrix) and returns
      a nX1 matrix called w.
    Example
      P := matrix{{.5, .25, .25}, {.5, 0, .5}, {.25, .25, .5}}
      findingw(P)

///

doc ///
  Key
    createW
    (createW, TransitionMatrix)
    (createW, Matrix)
  Usage
    createW(P)
  Description
    Text
      This function requires the input of P (the transition matrix) and returns
      a nXn matrix called W.
    Example
      P := matrix{{.5, .25, .25}, {.5, 0, .5}, {.25, .25, .5}}
      createW(P)

///

doc ///
  Key
    TransitionMatrix
    transitionMatrix
    (transitionMatrix, Matrix)
  Description
    Text
      The Transition Matrix is a square matrix where all of the rows sum up to 1.
    Example
      P := matrix{{.5, .25, .25}, {.5, 0, .5}, {.25, .25, .5}}
      transitionMatrix P
      Q := matrix{{.5, .25, .25}, {.5, 0, .5}}
      stopIfError = false
      transitionMatrix Q
      
///

doc ///
  Key
    ErgodicRegularMarkovChains
  Description
    Text
      This package helps the user determine if a Markov Chain is Regular, and what the transition matrix is approaching
///

-* Test section *-
TEST ///
  assert(isRegularMarkov (matrix{{.5, .5}, {.5, .5}},90));
  assert(not isRegularMarkov (matrix{{1,0,0},{0,1,0},{0,0,1}}, 5));
  
///

TEST ///
  assert zero clean(1e-15, findingw (matrix{{.5, .25, .25}, {.5, 0, .5}, {.25, .25, .5}}) - matrix{{.4, .2, .4}});
  assert not zero clean(1e-15, findingw (matrix{{.5, .25, .25}, {.5, 0, .5}, {.25, .25, .5}}) - matrix{{.5, .5, 0}});

///

TEST ///
  assert zero clean(1e-15, createW (matrix{{.5, .25, .25}, {.5, 0, .5}, {.25, .25, .5}}) - matrix{{.4, .2, .4}, {.4, .2, .4}, {.4, .2, .4}});
  assert not zero clean(1e-15, createW (matrix{{.5, .25, .25}, {.5, 0, .5}, {.25, .25, .5}}) - matrix{{.7, .3, 0}, {.7,.3, 0}, {.7, .3, 0}});
  
 ///
