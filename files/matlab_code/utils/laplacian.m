function L = laplacian(x)
    L = circshift(x,[1 0]) + circshift(x,[-1 0]) + circshift(x,[0 1]) + circshift(x,[0 -1]) - 4*x;
end
