class KIMCalculatorError(Exception):
  def __init__(self, msg):
    # Call the base class constructor with the parameters it needs
    super(KIMCalculatorError, self).__init__(msg)
    self.msg = msg
  def __str__(self):
    return self.msg
