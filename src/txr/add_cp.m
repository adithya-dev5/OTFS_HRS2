function xcp = add_cp(x, Mcp)
x = x(:);
if Mcp==0, xcp = x; return; end
L = length(x); assert(Mcp < L, 'Mcp must be less than signal length');
prefix = x(end-Mcp+1:end);
xcp = [prefix; x];
end
