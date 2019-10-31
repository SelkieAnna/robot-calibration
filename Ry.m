function Ry = Ry(q)
    Ry = [ cos(q), 0, sin(q), 0;
                0, 1,      0, 0;
          -sin(q), 0, cos(q), 0;
                0, 0,      0, 1];
end