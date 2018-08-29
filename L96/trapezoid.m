% Apply trapezoid rule for integration

function z = trapezoid(a, dx)

z = ( cumsum(a) - a./2).*dx;
z = z / max(z);
