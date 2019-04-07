from datetime import datetime

epoch = datetime(1950, 1, 1)
t = datetime(1956, 3, 2)
diff = t-epoch
print (diff.days * 24 * 3600 + diff.seconds)