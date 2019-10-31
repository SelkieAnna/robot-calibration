function Rz = Rz(q)
    Rz = [cos(q), -sin(q), 0, 0;
          sin(q),  cos(q), 0, 0;
               0,       0, 1, 0;
               0,       0, 0, 1];
end