def a(g):
  def b(x):
    return x + 1
  return b(g+20)

print a(25)