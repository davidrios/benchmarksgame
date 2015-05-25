from math import sqrt, pow
from os import paramStr
from strutils import parseInt

const PI = 3.141592653589793
const SOLAR_MASS = (4 * PI * PI)
const DAYS_PER_YEAR = 365.24

type
  Body = tuple
    name: string
    x, y, z, vx, vy, vz, mass: float64

var bodies = [
  ("Sun", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, SOLAR_MASS),
  ("Jupiter",
    4.84143144246472090e+00,
    -1.16032004402742839e+00,
    -1.03622044471123109e-01,
    1.66007664274403694e-03 * DAYS_PER_YEAR,
    7.69901118419740425e-03 * DAYS_PER_YEAR,
    -6.90460016972063023e-05 * DAYS_PER_YEAR,
    9.54791938424326609e-04 * SOLAR_MASS),
  ("Saturn",
    8.34336671824457987e+00,
    4.12479856412430479e+00,
    -4.03523417114321381e-01,
    -2.76742510726862411e-03 * DAYS_PER_YEAR,
    4.99852801234917238e-03 * DAYS_PER_YEAR,
    2.30417297573763929e-05 * DAYS_PER_YEAR,
    2.85885980666130812e-04 * SOLAR_MASS
  ),
  ("Uranus",
    1.28943695621391310e+01,
    -1.51111514016986312e+01,
    -2.23307578892655734e-01,
    2.96460137564761618e-03 * DAYS_PER_YEAR,
    2.37847173959480950e-03 * DAYS_PER_YEAR,
    -2.96589568540237556e-05 * DAYS_PER_YEAR,
    4.36624404335156298e-05 * SOLAR_MASS
  ),
  ("Neptune",
    1.53796971148509165e+01,
    -2.59193146099879641e+01,
    1.79258772950371181e-01,
    2.68067772490389322e-03 * DAYS_PER_YEAR,
    1.62824170038242295e-03 * DAYS_PER_YEAR,
    -9.51592254519715870e-05 * DAYS_PER_YEAR,
    5.15138902046611451e-05 * SOLAR_MASS
  )
]


proc printf(formatstr: cstring) {.header: "<stdio.h>", importc: "printf", varargs.}


proc advance(bodies: var openarray[Body], dt: float64) =
  echo("called")
  for i in countup(0, len(bodies) - 1):
    var b = bodies[i]
    for j in countup(i + 1, len(bodies) - 1):
      var b2 = bodies[j]
      let
        dx = b.x - b2.x
        dy = b.y - b2.y
        dz = b.z - b2.z
        distance = sqrt(dx * dx + dy * dy + dz * dz)
        mag = dt / (distance * distance * distance)

      b.vx -= dx * b2.mass * mag
      b.vy -= dy * b2.mass * mag
      b.vz -= dz * b2.mass * mag
      b2.vx += dx * b.mass * mag
      b2.vy += dy * b.mass * mag
      b2.vz += dz * b.mass * mag

      printf("%d, %d, %f, %f, %f, %f, %f\n", i, j, b.vx, b.vy, b.vz, b2.vx, b2.vy, b2.vz);

  for b in mitems(bodies):
    b.x += dt * b.vx;
    b.y += dt * b.vy;
    b.z += dt * b.vz;


proc energy(bodies: var openarray[Body]): float64 =
  for i in countup(0, len(bodies) - 1):
    let b = bodies[i]
    result += 0.5 * b.mass * (b.vx * b.vx + b.vy * b.vy + b.vz * b.vz)
    for j in countup(i + 1, len(bodies) - 1):
      let
        b2 = bodies[j]
        dx = b.x - b2.x
        dy = b.y - b2.y
        dz = b.z - b2.z
        distance = sqrt(dx * dx + dy * dy + dz * dz)
      result -= (b.mass * b2.mass) / distance


iterator mpairsp1[T](a: openarray[T]): tuple[a1, a2: T] =
  for i in countup(0, len(a) - 1):
    for j in countup(i + 1, len(a) - 1):
      yield (a[i], a[j])


proc energy2(bodies: var openarray[Body]): float64 =
  result = 0.0
  for b, b2 in mpairsp1(bodies):
    let
      dx = b.x - b2.x
      dy = b.y - b2.y
      dz = b.z - b2.z
    result -= ((b.mass * b2.mass) / pow(dx * dx + dy * dy + dz * dz, 0.5))

  for b in items(bodies):
    result += (b.mass * (b.vx * b.vx + b.vy * b.vy + b.vz * b.vz) / 2.0)


proc offset_momentum(bodies: var openarray[Body]) =
  var
    px = 0.0
    py = 0.0
    pz = 0.0

  for body in items(bodies):
    px += body.vx * body.mass;
    py += body.vy * body.mass;
    pz += body.vz * body.mass;

  bodies[0].vx = - px / SOLAR_MASS;
  bodies[0].vy = - py / SOLAR_MASS;
  bodies[0].vz = - pz / SOLAR_MASS;


when isMainModule:
  let n = parseInt(1.paramStr.string)

  offset_momentum(bodies)
  printf("%.9f\n", energy(bodies))

  for i in 0..n-1:
    advance(bodies, 0.01)

  printf("%.9f\n", energy(bodies))