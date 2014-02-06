function y = apply_deflprec(x, M1, T, U)

y = M1*(x + U*(T\(U'*x)));

end
