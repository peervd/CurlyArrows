# Logging Configuratie

De CurlyArrows applicatie heeft een uitgebreid logging systeem dat zowel naar console als naar bestanden schrijft.

## Overzicht

Het logging systeem biedt:
- **Console logging** - Altijd ingeschakeld voor realtime monitoring
- **File logging** - Optioneel, naar configureerbare directory
- **Automatische rotatie** - Logbestanden worden automatisch geroteerd bij 10MB
- **Verschillende log levels** - DEBUG, INFO, WARNING, ERROR, CRITICAL
- **Aparte error logs** - Fouten worden ook apart opgeslagen

## Configuratie

### Omgevingsvariabelen

```bash
# Log directory (default: /tmp/logs voor productie, ./logs voor development)
LOG_DIR=/tmp/logs

# Log level (default: INFO voor productie, DEBUG voor development)
LOG_LEVEL=INFO
```

### Configuratie per environment

#### Development
```bash
export FLASK_ENV=development
export LOG_DIR=./logs
export LOG_LEVEL=DEBUG
```

#### Production (Azure Container Apps)
```bash
export FLASK_ENV=production
export LOG_DIR=/tmp/logs
export LOG_LEVEL=INFO
```

## Logbestanden

Wanneer file logging is ingeschakeld, worden twee bestanden aangemaakt:

### 1. curlyarrows.log
Bevat alle log berichten vanaf het geconfigureerde level (standaard INFO):
- Applicatie startup informatie
- User acties
- API calls
- Database operaties
- Algemene informatie

### 2. curlyarrows_errors.log
Bevat alleen ERROR en CRITICAL berichten:
- Exception tracebacks
- Database errors
- Authentication failures
- API errors

## Log Rotatie

Logbestanden worden automatisch geroteerd:
- **Maximum grootte**: 10MB per bestand
- **Backup count**: 5 bestanden
- **Naamgeving**: `curlyarrows.log.1`, `curlyarrows.log.2`, etc.

## Gebruik in Code

### In Route Handlers (Blueprints)

```python
from flask import current_app

@app.route('/example')
def example():
    current_app.logger.info("User accessed example page")

    try:
        # Your code here
        result = do_something()
        current_app.logger.debug(f"Result: {result}")
    except Exception as e:
        current_app.logger.error(f"Error in example: {e}", exc_info=True)

    return "OK"
```

### In Andere Modules

```python
from logging_config import get_logger

logger = get_logger(__name__)

def process_data(data):
    logger.info(f"Processing data: {len(data)} items")

    try:
        # Process data
        logger.debug(f"Data details: {data}")
    except Exception as e:
        logger.error(f"Failed to process data: {e}", exc_info=True)
        raise
```

## Log Levels

### DEBUG (10)
Gedetailleerde informatie voor debugging:
```python
app.logger.debug("Variable x = " + str(x))
```

### INFO (20)
Algemene informatie over applicatie werking:
```python
app.logger.info("User logged in successfully")
```

### WARNING (30)
Waarschuwingen die aandacht vereisen maar geen direct probleem:
```python
app.logger.warning("API rate limit approaching")
```

### ERROR (40)
Fouten die de werking beïnvloeden:
```python
app.logger.error("Database connection failed", exc_info=True)
```

### CRITICAL (50)
Kritieke fouten die de applicatie kunnen stoppen:
```python
app.logger.critical("Configuration file missing!")
```

## Azure Container Apps

### Logs Bekijken

```bash
# Recente logs (console output)
az containerapp logs show --name curlyarrows-app --resource-group eur-rg-curly --tail 100

# Logs realtime volgen
az containerapp logs show --name curlyarrows-app --resource-group eur-rg-curly --follow

# Alleen errors
az containerapp logs show --name curlyarrows-app --resource-group eur-rg-curly --tail 50 | grep ERROR
```

### Log Persistentie

**Belangrijk**: Logs in `/tmp/logs` zijn tijdelijk en gaan verloren bij:
- Container herstart
- Nieuwe deployment
- Scale down naar 0 replicas

### Oplossingen voor Persistente Logs

#### Optie 1: Azure Application Insights (Aanbevolen)
```python
# Voeg toe aan requirements.txt
opencensus-ext-azure
opencensus-ext-flask

# In app.py
from opencensus.ext.azure.log_exporter import AzureLogHandler
handler = AzureLogHandler(connection_string='InstrumentationKey=...')
app.logger.addHandler(handler)
```

#### Optie 2: Azure Storage (File Share)
Mount een Azure File Share voor persistente logs:
```bash
az containerapp create \
  --storage-account myaccount \
  --storage-name logs \
  --storage-type AzureFile
```

#### Optie 3: Centralized Logging Service
- Elastic Stack (ELK)
- Azure Monitor
- Datadog
- Splunk

## Testen

Test de logging configuratie lokaal:

```bash
# Maak logs directory
mkdir -p logs

# Run test script
python test_logging.py

# Bekijk de logs
cat logs/curlyarrows.log
cat logs/curlyarrows_errors.log
```

## Best Practices

### 1. Gebruik het juiste log level
```python
# ✅ GOED
app.logger.info("User logged in")
app.logger.error("Database query failed", exc_info=True)

# ❌ FOUT
app.logger.error("User logged in")  # Dit is geen error
app.logger.info("Critical system failure")  # Dit zou CRITICAL moeten zijn
```

### 2. Voeg context toe
```python
# ✅ GOED
app.logger.error(f"Failed to process order {order_id} for user {user_id}", exc_info=True)

# ❌ FOUT
app.logger.error("Failed to process order")  # Te weinig context
```

### 3. Gebruik exc_info voor exceptions
```python
# ✅ GOED
try:
    dangerous_operation()
except Exception as e:
    app.logger.error(f"Operation failed: {e}", exc_info=True)

# ❌ FOUT
except Exception as e:
    app.logger.error(str(e))  # Mist stack trace
```

### 4. Vermijd gevoelige informatie
```python
# ✅ GOED
app.logger.info(f"User {user_id} logged in")

# ❌ FOUT - bevat wachtwoord!
app.logger.debug(f"Login attempt: {username} / {password}")
```

### 5. Log bij belangrijke events
- User authentication (login/logout)
- Database operaties (vooral writes)
- API calls (intern en extern)
- Error conditions
- Security events

## Troubleshooting

### Logs worden niet geschreven naar bestand

**Probleem**: Alleen console output, geen logbestanden

**Oplossing**:
```bash
# Controleer of LOG_DIR bestaat en schrijfbaar is
mkdir -p $LOG_DIR
chmod 755 $LOG_DIR

# Test of Python kan schrijven
python -c "import os; print(os.access('$LOG_DIR', os.W_OK))"
```

### Te veel log output

**Probleem**: Logbestanden worden te groot

**Oplossing**:
```bash
# Verhoog log level naar WARNING of ERROR
export LOG_LEVEL=WARNING

# Of in config.py:
LOG_LEVEL = 'WARNING'
```

### Geen logs in Azure

**Probleem**: `az containerapp logs show` toont niets

**Oplossing**:
```bash
# Controleer of container draait
az containerapp show --name curlyarrows-app --resource-group eur-rg-curly --query properties.runningStatus

# Controleer recente events
az containerapp logs show --name curlyarrows-app --resource-group eur-rg-curly --type system
```

## Voorbeeld Output

### Console (Development)
```
[12:34:56] INFO: ============================================================
[12:34:56] INFO: CurlyArrows Application Starting
[12:34:56] INFO: Environment: development
[12:34:56] INFO: Log Level: DEBUG
[12:34:56] INFO: Debug Mode: True
[12:34:56] INFO: Log Directory: ./logs
[12:34:56] INFO: ============================================================
[12:34:56] INFO: Creating application with config: development
[12:34:56] INFO: Flask-Session initialized
[12:34:56] INFO: Blueprints registered: main, auth, api
[12:34:56] INFO: Database initialized
```

### File (curlyarrows.log)
```
[2025-12-03 12:34:56] INFO in app (create_app:23): Creating application with config: development
[2025-12-03 12:34:56] INFO in app (create_app:27): Flask-Session initialized
[2025-12-03 12:34:56] INFO in app (create_app:37): Blueprints registered: main, auth, api
[2025-12-03 12:34:56] INFO in app (create_app:43): Database initialized
[2025-12-03 12:35:10] INFO in main (index:16): User John Doe accessed index page
```

## Meer Informatie

- [Python Logging Documentation](https://docs.python.org/3/library/logging.html)
- [Flask Logging](https://flask.palletsprojects.com/en/2.3.x/logging/)
- [Azure Container Apps Logging](https://learn.microsoft.com/en-us/azure/container-apps/logging)
