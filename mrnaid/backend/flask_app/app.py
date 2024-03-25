import sys
import os

# Add the directory containing Evaluation.py to the PYTHONPATH
sys.path.append(os.path.abspath("../common"))
print(sys.executable)
from .routes import app

if __name__ == "__main__":
    app.run(debug=True)