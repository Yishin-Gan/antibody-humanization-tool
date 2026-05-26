# Workstation Runbook — Running the Antibody Humanization Advisor

This is the **everyday operating manual** for the person who inherits this
specific workstation. You don't need to install anything — everything is
already set up. You just need to start it correctly each morning and know
what to do when something looks wrong.

If you've never used a terminal before, **read this whole document once
before you try anything.** It's all copy-paste — you don't need to type
or improvise.

---

## At-a-glance cheat sheet

For days when you just need to start the app:

1. From the Windows Start menu, search for **"Ubuntu"** and open it. A black-on-white terminal window appears.
2. Copy and paste the following block exactly (right-click in the terminal to paste):
   ```bash
   cd /home/yishin/humanization
   docker ps | grep -q humanization || ./docker_run.sh
   sleep 3
   CONTAINER=$(docker ps --format '{{.Names}}' | head -1)
   docker exec -it "$CONTAINER" bash -c "cd /workspace/antibody-humanization-tool && python3 run_web.py --host 0.0.0.0 --port 8080"
   ```
3. Wait until you see lines that include `* Running on http://0.0.0.0:8080`. The app is now live.
4. From any browser on the office network: open `http://10.6.32.111:8080`
5. **Do not close the Ubuntu window.** Minimize it. Closing it shuts down the app.

That's the whole daily operation. Everything below is detail, troubleshooting,
and "what to do if step 3 doesn't appear".

---

## Table of contents

1. [Background — what you are starting](#1-background--what-you-are-starting)
2. [First-time setup — do this once](#2-first-time-setup--do-this-once)
3. [Starting the app each morning — detailed walk-through](#3-starting-the-app-each-morning--detailed-walk-through)
4. [Accessing the app from any browser](#4-accessing-the-app-from-any-browser)
5. [Sharing the URL with the team](#5-sharing-the-url-with-the-team)
6. [Stopping the app at the end of the day](#6-stopping-the-app-at-the-end-of-the-day)
7. [Common problems and fixes](#7-common-problems-and-fixes)
8. [What NOT to do](#8-what-not-to-do)
9. [If everything is broken — start over from scratch](#9-if-everything-is-broken--start-over-from-scratch)
10. [Asking for help](#10-asking-for-help)

---

## 1. Background — what you are starting

The Antibody Humanization Advisor is a web app. The way it runs on this
workstation is a sandwich of layers:

```
┌─────────────────────────────────────────────────┐
│  Windows 11 (this workstation)                  │
│  ┌───────────────────────────────────────────┐  │
│  │  WSL — Ubuntu 22.04                       │  │
│  │  (this is where you type commands)        │  │
│  │  ┌─────────────────────────────────────┐  │  │
│  │  │  Docker container (nvidia/pytorch)  │  │  │
│  │  │  ┌──────────────────────────────┐   │  │  │
│  │  │  │  python3 run_web.py          │   │  │  │
│  │  │  │  ← the actual app            │   │  │  │
│  │  │  └──────────────────────────────┘   │  │  │
│  │  └─────────────────────────────────────┘  │  │
│  └───────────────────────────────────────────┘  │
└─────────────────────────────────────────────────┘
```

To start the app, you have to go inward through the layers. You open
Ubuntu, then start the container, then start the Python program inside
the container. The cheat sheet above does all three with copy-paste.

**You only ever need to touch the Ubuntu terminal.** Never go into
Docker Desktop directly, never open PowerShell, never click into VS Code.

---

## 2. First-time setup — do this once

You should not need to do anything here on day one — the prior maintainer
left everything installed and configured. This section is so you can
confirm everything is in place before you trust the workstation for daily
use.

### 2.1 Log in to the workstation

Use the Windows account you've been told to use for this role. **It should
be a shared "lab" or service account, not your personal account** — that
way the setup survives staff turnover. If you were given your own personal
Windows account for this, ask IT to switch you to the shared account
before you go further. (Reason: when your personal account is eventually
deleted, Docker and all its containers get deleted with it.)

### 2.2 Confirm Ubuntu is installed

1. Click the Windows Start button (bottom-left corner).
2. Type "Ubuntu" and look for an entry called **Ubuntu** or **Ubuntu 22.04** (the icon is an orange circle).
3. Click it. A black-on-white terminal window opens.
4. You should see a prompt like `yishin@R407-0003196:~$` or similar (your username will differ).

If Ubuntu does not appear in the Start menu, ask IT to install WSL with
Ubuntu 22.04. The successor cannot run this app without it.

### 2.3 Confirm Docker Desktop is running

1. Look at the system tray (bottom-right corner of the Windows screen, near the clock). You may need to click the small `^` to see hidden icons.
2. Look for a small **blue whale icon**. Hover over it.
3. The tooltip should say "Docker Desktop is running".

If the whale is grey or not present:
- Click the Windows Start button, type "Docker Desktop", and click it.
- Wait 30-60 seconds. The whale icon appears and eventually turns blue.
- If Docker Desktop is not installed at all, ask IT. The successor cannot run this app without it.

### 2.4 Confirm the project files are in place

In the Ubuntu terminal you opened in step 2.2, type:

```bash
ls /home/yishin/humanization/
```

You should see at minimum these items:

```
antibody-humanization-tool   data   docker_run.sh
```

If anything is missing (especially `docker_run.sh` or `data/`), **stop
and ask for help before continuing.** Re-creating these from scratch is
a 1-hour install, not a daily operation.

Also verify the OASis humanness database is present:

```bash
ls -lh /home/yishin/humanization/data/OASis_9mers_v1.db
```

You should see a file size of about **23G** (23 gigabytes). If the file
is missing, the app will start but every submission will fail. Re-downloading
it takes ~30 minutes; see Section 9.

---

## 3. Starting the app each morning — detailed walk-through

This is the "long version" of the cheat sheet at the top, with what to
expect at each step.

### 3.1 Open the Ubuntu terminal

Same as Section 2.2. The terminal window has a black-on-white background
and shows a prompt like `yishin@R407-0003196:~$`.

### 3.2 Check if the Docker container is already running

It might already be running from a previous day. Type:

```bash
docker ps
```

Press Enter. You'll see one of two things:

- **Table with at least one row** that mentions `nvidia/pytorch` or has a NAMES column with a value → container is already running, skip to step 3.4.
- **Header row only, no data rows** → container is not running, continue to step 3.3.

### 3.3 Start the Docker container

Type:

```bash
cd /home/yishin/humanization
./docker_run.sh
```

You'll see a few seconds of Docker output, ending with a long alphanumeric
string (the container ID). That means it started in the background.

Confirm it's running:

```bash
docker ps
```

A row should now appear in the table.

### 3.4 Get the container's name

```bash
docker ps --format '{{.Names}}'
```

The output is one or more words. The first word is your container name —
you'll use it in the next step. Call this value `CONTAINER_NAME`.

For example, if the output is:

```
pytorch_dev
```

then `CONTAINER_NAME` is `pytorch_dev`.

### 3.5 Enter the container and start the Flask app

Type the following, replacing `<CONTAINER_NAME>` with the actual name
from step 3.4:

```bash
docker exec -it <CONTAINER_NAME> bash
```

The prompt changes — usually to something like `root@docker-desktop:/workspace#`.
You are now inside the container.

Navigate to the project folder and start the app:

```bash
cd /workspace/antibody-humanization-tool
python3 run_web.py --host 0.0.0.0 --port 8080
```

After about 5 seconds, you should see lines that include:

```
Serving Antibody Humanization Advisor on http://0.0.0.0:8080
 * Running on http://0.0.0.0:8080
 * Running on http://127.0.0.1:8080
Press CTRL+C to quit
```

**The app is now running.** The terminal will stay paused on this output —
that is normal. The Flask app is running and waiting for browser requests.

### 3.6 IMPORTANT — leave the terminal window open

Do **not** close this Ubuntu window. Do **not** press Ctrl+C. Do **not**
press Ctrl+D. Any of those will shut down the app.

You can **minimize** the window (the `_` button in the top-right). Just
don't close it.

---

## 4. Accessing the app from any browser

The app is now reachable at two URLs:

| From where | URL to use |
|---|---|
| **The same workstation** (the one running the app) | `http://localhost:8080` |
| **Any other computer on the office network** | `http://10.6.32.111:8080` |

To use the app, open any web browser (Chrome, Edge, Firefox, Safari) and
type the URL into the address bar. The input form should load.

If the workstation's IP changes (rare — usually only after IT reconfigures
the network), the second URL changes too. To find the current IP at any
time, run in PowerShell on the workstation:

```powershell
ipconfig | findstr IPv4
```

Look for the entry under "Ethernet adapter" (not Wi-Fi, not WSL, not
Bluetooth). That number is the new IP.

---

## 5. Sharing the URL with the team

When you want others to use the tool, send them this template message:

> The antibody humanization tool is at:
>
> http://10.6.32.111:8080
>
> Open it in any browser (Chrome, Edge, Firefox). You need to be on the
> office network — VPN-in from home works, public Wi-Fi does not.
>
> One submission takes 1–2 minutes. Concurrent submissions queue, so if
> two of you submit at once, the second waits for the first to finish.

---

## 6. Stopping the app at the end of the day

You have three options. Pick the one that matches your situation.

### Option A — Leave it running overnight (recommended)

Do nothing. Lock the Windows screen (Win+L), but don't close the Ubuntu
window or shut down the workstation. The app stays available to the team
for early-morning users.

Resource usage when idle: roughly 2 GB of RAM, near-zero CPU. The
workstation will not slow down.

### Option B — Stop the app but leave the workstation on

In the Ubuntu terminal where the app is running, press **Ctrl+C** once.
You'll see the app print a shutdown message and the prompt returns:

```
^C
root@docker-desktop:/workspace/antibody-humanization-tool#
```

The container is still running, just not serving the web app. To restart
the app later without restarting the container, just re-run:

```bash
python3 run_web.py --host 0.0.0.0 --port 8080
```

### Option C — Full shutdown (only if you're going away for a week+)

Stop the app (Ctrl+C as above), exit the container (`exit` then Enter),
then stop the container:

```bash
docker stop $(docker ps --format '{{.Names}}' | head -1)
```

After this, the next morning you have to do the full startup from
Section 3.

---

## 7. Common problems and fixes

### "I see 'permission denied' or 'cannot connect to docker daemon'"

Docker Desktop isn't running. Click the Windows Start button, type
"Docker Desktop", click it, wait 30 seconds for the whale icon to turn
blue, then retry your command.

### "I closed the Ubuntu window by accident"

The app is dead. Start over from Section 3.

### "The browser says ERR_CONNECTION_REFUSED"

The most likely cause: the Flask app isn't actually running. Check the
Ubuntu terminal where you started it. If it shows a normal `$` prompt
instead of the `* Running on http://...` lines, the app stopped. Restart
from Section 3.5.

Second most likely cause: the workstation's IP changed. Re-check the IP
per Section 4 and use the new one.

### "The page loads but submissions hang or fail with 'OASis DB not found'"

The OASis humanness database is missing from the expected location.
Check:

```bash
ls -lh /home/yishin/humanization/data/OASis_9mers_v1.db
```

Should show a ~23 GB file. If it says "No such file or directory", the
file is missing and needs to be restored. See Section 9.

### "The first submission after starting takes forever (30+ seconds)"

That's normal. The Sapiens model weights load into memory on the first
request. Subsequent requests are fast. You can pre-warm by submitting
one short sequence before the team is watching.

### "The team says the URL doesn't work from their laptops"

Check three things in order:
1. Is the workstation on? (you're using it, so yes)
2. Is the app running? (the Ubuntu terminal shows `* Running on http://...`)
3. Is everyone on the same office network? Wi-Fi guests cannot reach
   the LAN.

### "After a Windows update / reboot, nothing works"

Expected. A reboot wipes the running container. Just do Section 3 from
the start — `./docker_run.sh` will create a fresh container, the app
starts inside it, life continues.

### "I see scary red text I don't understand"

Don't panic. Take a screenshot, then send it to whoever is your
technical contact (Section 10). Most Python tracebacks look terrifying
but are routine. Do not type random commands trying to "fix" things —
that's how an already-bad day becomes a much worse one.

---

## 8. What NOT to do

Things that will break the setup and require asking for help:

- **Don't run `rm -rf` or `delete`** anything inside `/home/yishin/humanization/`. The OASis database alone took 30 minutes to download and the source code is the entire project.
- **Don't `git pull`** unless you've been told to. Code updates can introduce changes you haven't tested.
- **Don't `pip install`** anything inside the container. The Python environment was carefully assembled. New packages can break existing ones in non-obvious ways.
- **Don't reboot the workstation** unless Windows Update forces it. The container goes away on reboot and you have to restart it (Section 3).
- **Don't shut down Docker Desktop** unless you mean to. Stopping Docker stops the container, which stops the app.
- **Don't change the workstation's IP** by hand or change Wi-Fi/Ethernet settings. The team's URL bookmarks will break.
- **Don't share the URL with anyone outside the company.** There is no authentication on the app — anyone with the URL on the network can use it.

---

## 9. If everything is broken — start over from scratch

If the workstation is in such a bad state that the steps in Section 3 no
longer work, you have a nuclear option: re-install the tool from the
GitHub repository. This takes about 1 hour and gives you a fresh,
guaranteed-working install.

```bash
# Open Ubuntu terminal
cd /home/yishin/humanization
# (Optional but recommended: back up the existing folder first)
mv antibody-humanization-tool antibody-humanization-tool.broken.$(date +%Y%m%d)
# Get a fresh copy of the code
git clone https://github.com/Yishin-Gan/antibody-humanization-tool.git
cd antibody-humanization-tool
# Make sure the OASis database is in the right place
ln -sf /home/yishin/humanization/data/OASis_9mers_v1.db data/OASis_9mers_v1.db
# Run the proper installer (this builds a clean container)
./install.sh
```

After `install.sh` finishes, the tool runs on port **5000** (not 8080)
via the proper production container. URL becomes `http://10.6.32.111:5000`.
Tell the team about the new port.

If `install.sh` fails, you cannot fix it yourself — escalate (Section 10).

---

## 10. Asking for help

When something is wrong, the more information you provide, the faster
it gets fixed. Send a message with:

1. **What you were trying to do** — "I wanted to start the app this morning."
2. **What happened instead** — "The browser shows 'site can't be reached'."
3. **Screenshot of the Ubuntu terminal** showing the last 20 lines.
4. **Screenshot of the browser error.**
5. **Output of these diagnostic commands**, copy-pasted from Ubuntu:

   ```bash
   docker ps
   docker info | grep "Server Version"
   ls /home/yishin/humanization/
   ls /home/yishin/humanization/data/OASis_9mers_v1.db
   curl http://localhost:8080
   ```

Send to: **[fill in technical contact name and email/Slack here]**

Or open a GitHub issue at:
https://github.com/Yishin-Gan/antibody-humanization-tool/issues

---

## Quick-reference card (print this and tape it to the workstation)

```
START THE APP:
  1. Start menu → Ubuntu
  2. cd /home/yishin/humanization
  3. ./docker_run.sh         (skip if already running — check with "docker ps")
  4. docker exec -it <NAME-FROM-docker-ps> bash
  5. cd /workspace/antibody-humanization-tool
  6. python3 run_web.py --host 0.0.0.0 --port 8080
  7. Wait for "* Running on http://0.0.0.0:8080"
  8. Browser: http://10.6.32.111:8080
  9. LEAVE THE WINDOW OPEN.

STOP THE APP:
  In the Ubuntu window: press Ctrl+C once.

SOMETHING IS BROKEN:
  See WORKSTATION_RUNBOOK.md Section 7 (Common problems).
  Or contact: [name/email]
```
