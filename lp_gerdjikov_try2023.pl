%set_prolog_flag(answer_write_options,[max_depth(0)]).
nat(c).
nat(s(X)):-nat(X).
%X+0=X
sum(X,c,X).
%X+(Y+1)=(X+Y)+1
sum(X,s(Y),s(Z)):-sum(X,Y,Z).
%X<Y <-> exists Z\in N (X+(Z+1)=Y)
less(X,Y):-sum(X,s(Z),Y).
less_eq(X,Y):-sum(X,_,Y).
%X.0=0; X.(Y+1)=(X.Y)+X
prod(X,c,c).
prod(X,s(Y),Z):-prod(X,Y,Z1),sum(Z1,X,Z).
%X|Z+1 <-> X\in {0,1,....,Z+1} & \exists Y\in {0,1,...,Z+1} (X.Y=Z+1)
divisors(X,s(Z)):-less_eq(X,s(Z)),less_eq(Y,s(Z)),prod(X,Y,s(Z)).
divisors(X,c):-nat(X).

%Bezout: ua+bv=gcd(a,b)
%(a,b)|->(u,v,d)
bezout(0,B,0,1,B).
%(b-a,b,u1,v,d): au1+(b-a)v=d, then a(u1-v) + bv=d, i.e. u=u1-v.
bezout(A,B,U,V,D):-B>=A,A>0,B1 is B-A,bezout(A,B1,U1,V,D),U is U1 - V.
% N x N = {0} x N \cup {(A,B)\in N x N | 0<A<=B} \cup {(A,B) \in N x N |
% B<A}
bezout(A,B,U,V,D):-A>B,bezout(B,A,V,U,D).

fibonacci(N,N):-N=<1.
fibonacci(F,N):-N>1,N1 is N-1,fibonacci(F1,N1),N2 is N-2,fibonacci(F2,N2),F is F1 + F2.
%N -> FN,F{N+1}
fibonacci1(0,0,1).
fibonacci1(N,F,F_next):-N>0,N1 is N-1,fibonacci1(N1,F_prev,F),F_next is F + F_prev.
matrix_mult(A,B,C,D,A1,B1,C1,D1,A2,B2,C2,D2):-A2 is A*A1 + B*C1,B2 is A*B1+B*D1,C2 is C*A1+D*C1,D2 is C*B1+D*D1.
matrix_fast(0,A,B,C,D,1,0,0,1).
matrix_fast(1,A,B,C,D,A,B,C,D).
matrix_fast(N,A,B,C,D,AN,BN,CN,DN):-N>1,N mod 2=:=0,N1 is N/2,matrix_fast(N1,A,B,C,D,A1,B1,C1,D1),matrix_mult(A1,B1,C1,D1,A1,B1,C1,D1,AN,BN,CN,DN).
matrix_fast(N,A,B,C,D,AN,BN,CN,DN):-N>1,N mod 2=\=0,N1 is N-1,matrix_fast(N1,A,B,C,D,A1,B1,C1,D1),matrix_mult(A,B,C,D,A1,B1,C1,D1,AN,BN,CN,DN).
fast_fibonacci(N,F):-matrix_fast(N,0,1,1,1,_,F,_,_).

%list L=[H|T], T - list; L=[].
%#=
len([],0).
len([H|T],N):-len(T,N1),N is N1+1.
is_eq_len([],0).
is_eq_len([H|T],N):-N1 is N-1,is_eq_len(T,N1).
eq_len(L,N):-len(L,N1),N1=:=N.

split([],[],[]).
split([H],[H],[]).
split([A,B|T],[A|T_even],[B|T_odd]):-split(T,T_even,T_odd).

merge1([],SL2,SL2).
merge1([H1|T1],[],[H1|T1]).
merge1([H1|T1],[H2|T2],[H1|T]):-H1=<H2,merge1(T1,[H2|T2],T).
merge1([H1|T1],[H2|T2],[H2|T]):-H2<H1,merge1([H1|T1],T2,T).

merge_sort([],[]).
merge_sort([H],[H]).
merge_sort([A,B|T],SL):-split([A,B|T],L1,L2),merge_sort(L1,SL1),merge_sort(L2,SL2),merge1(SL1,SL2,SL).

%given L1 and L2 generates L=L1.L2
%and conversely, given L, generates all (L1,L2) s.t. L1.L2=L
concat1([],L,L).
concat1([H|T],L2,[H|T2]):-concat1(T,L2,T2).
% is_sorted(L) checks whether a list of integers is sorted in increasing
% order
is_not_sorted(L):-concat1(_,[A,B|_],L),A>B.
is_sorted(L):-not(is_not_sorted(L)).
% mem2(L,X), given a list L, generates (under oversatisfaction) in X the
% elements of L.
mem2(L,X):-concat1(_,[X|_],L).

perm([],[]).
perm([A|T],P):-perm(T,PT),concat1(P1,P2,PT),concat1(P1,[A|P2],P).

my_sort(L,S):-perm(L,S),is_sorted(S).

subseq([],[]).
subseq([A|T],[A|S]):-subseq(T,S).
subseq([_|T],S):-subseq(T,S).

subseq1(_,[]).
subseq1(L,[A|S]):-concat1(_,[A|T],L),subseq1(T,S).

% gen_incr_subseq(L,IS) - given L, generates in IS all increasing
% subsequences of L
gen_incr_subseq(L,IS):-subseq(L,IS),is_sorted(IS).
% not_max_incr_subseq(L,N) checks whethere there is an increasing
% subsequence of L with length greater than N
not_max_incr_subseq(L,N):-gen_incr_subseq(L,IS),len(IS,N1),N1>N.
max_incr_subseq(L,IS):-gen_incr_subseq(L,IS),len(IS,N),not(not_max_incr_subseq(L,N)).

is_edge(U,V,E):-mem2(E,[U,V]).
is_edge(U,V,E):-mem2(E,[V,U]).




























