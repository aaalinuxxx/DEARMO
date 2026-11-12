FROM python:3.11-slim

WORKDIR /app

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/app

RUN apt-get update && apt-get install -y --no-install-recommends \
    libxrender1 \
    libxext6 \
    libsm6 \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY . /app

RUN mkdir -p /app/outputs/exports /app/outputs/demo /app/data/cache/pubchem /app/data/cache/chembl /app/data/cache/bindingdb /app/data/cache/pdb

ENTRYPOINT ["python", "cli.py"]
CMD ["demo", "--export-dir", "outputs/demo"]
