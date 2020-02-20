from calculateS import calculateSCoordinate

h = 5
print(calculateSCoordinate(0.07833071,0.09612648))

z = [ 0.124,       0.07833071, -0.0240309,  -0.12887921, -0.23372753, -0.33857584,-0.33857584, -0.23372753, -0.12887921, -0.0240309,   0.07833071]
y = [ 0.,          0.09612648,  0.11637895,  0.08312782,  0.04987669,  0.01662556,-0.01662556, -0.04987669, -0.08312782, -0.11637895, -0.09612648]

sList = []
nList = []

for i in range(len(z)):
    s,n = calculateSCoordinate(z[i],y[i])
    sList.append(s)
    nList.append(n)


print(sList)
print(nList)
