# Version 1.0.1 - 2026-04-27

The previously fixed batch size for logging ANIm, dnadiff and ANIb comparisons
every 50 comparisons was made dynamic to log more often with larger genomes.
This is important if the worker nodes get interrupted without being able to save
completed but unrecorded results.

# Version 1.0.0 - 2025-08-22

The initial public release of pyANI-plus.
