# Computer Networks (CN): Comprehensive Study Guide
## All Numerical Problems with Solutions, Theoretical Q&A / Differences, and Code Implementations
**Course**: 25AID202 - Introduction to Computer Networks (Kurose & Ross 8th/9th Edition)  
**Coverage**: Chapters 1 through 8 (All PPT Slides)

---

## TABLE OF CONTENTS
1. [PART 1: ALL NUMERICAL PROBLEMS & STEP-BY-STEP SOLUTIONS](#part-1-all-numerical-problems--step-by-step-solutions)
   - [Chapter 1: Transmission Delay, Propagation Delay, Queuing Delay, End-to-End Delay, Packet vs. Circuit Switching, Throughput, Topology Formulas](#chapter-1-numericals)
   - [Chapter 2: HTTP Response Times, Web Caching Link Utilization & Delay, P2P vs. Client-Server Distribution Times](#chapter-2-numericals)
   - [Chapter 3: UDP Checksum, Stop-and-Wait Utilization, Pipelining Efficiency, Selective Repeat Window Sizes, TCP RTT & Timeout Estimation, TCP Throughput & Long-Fat Pipes](#chapter-3-numericals)
   - [Chapter 4: Router Buffering Sizing, Classful IP Addressing, CIDR Subnetting & Address Allocations, NAT Translations & Port Mappings, IPv4 Fragmentation & Reassembly](#chapter-4-numericals)
   - [Chapter 5: Dijkstra's Link-State Algorithm Step-by-Step, Bellman-Ford Distance Vector Matrix Iterations, Count-to-Infinity Escalation](#chapter-5-numericals)
   - [Chapter 6: Parity Bits (1D & 2D), CRC (Cyclic Redundancy Check) Division, Slotted ALOHA & Pure ALOHA Efficiency, CSMA/CD Efficiency, Switch Self-Learning Tables](#chapter-6-numericals)
   - [Chapter 7: SNR (dB) Calculations, Shannon Channel Capacity, Digital Modulation Bits/Symbol, 4G/5G Resource Blocks](#chapter-7-numericals)
   - [Chapter 8: Modular Arithmetic, RSA Public-Key Cryptography (Key Generation, Encryption, Decryption), Internet Checksum Hash Collisions](#chapter-8-numericals)
2. [PART 2: THEORETICAL QUESTIONS, "WHY", AND COMPARISONS ("DIFFERENCE BETWEEN")](#part-2-theoretical-questions-why-and-comparisons)
   - [Chapter 1: Foundations, Protocol Stacks, Switching & Topologies](#chapter-1-theory)
   - [Chapter 2: Application Layer, Protocols & Architectures](#chapter-2-theory)
   - [Chapter 3: Transport Layer, Reliability & Congestion Control](#chapter-3-theory)
   - [Chapter 4: Network Layer (Data Plane), Addressing & Forwarding](#chapter-4-theory)
   - [Chapter 5: Network Layer (Control Plane), Routing & SDN](#chapter-5-theory)
   - [Chapter 6: Link Layer, Error Control, LANs & Multiple Access](#chapter-6-theory)
   - [Chapter 7: Wireless, Mobile Networks & Cellular/Satellite](#chapter-7-theory)
   - [Chapter 8: Network Security, Cryptography & Authentication](#chapter-8-theory)
3. [PART 3: ALL CODE SNIPPETS & ALGORITHMIC PSEUDOCODE](#part-3-all-code-snippets--algorithmic-pseudocode)
   - [Python UDP Client & UDP Server](#code-udp)
   - [Python TCP Client & TCP Server](#code-tcp)
   - [Socket Timeouts & Non-blocking I/O](#code-socket-events)
   - [Dijkstra's Link-State Algorithm Implementation](#code-dijkstra)
   - [Bellman-Ford Distance Vector Algorithm Implementation](#code-bellman-ford)
   - [CSMA/CD Binary Exponential Backoff Algorithm](#code-csmacd)

---

# PART 1: ALL NUMERICAL PROBLEMS & STEP-BY-STEP SOLUTIONS

<a name="chapter-1-numericals"></a>
## CHAPTER 1: INTRODUCTION & NETWORK PERFORMANCE

### Formula Summary
$$	ext{Transmission Delay } (d_{trans}) = \frac{L}{R} \quad (L = \text{packet length in bits}, R = \text{transmission rate in bps})$$
$$	ext{Propagation Delay } (d_{prop}) = \frac{d}{s} \quad (d = \text{link distance in meters}, s = \text{propagation speed } \approx 2 \times 10^8 \text{ m/s in fiber/copper})$$
$$	ext{Nodal Processing Delay } (d_{proc}) = \text{time to check bit errors and determine output link}$$
$$	ext{Queueing Delay } (d_{queue}) = \text{time waiting at output buffer for transmission}$$
$$d_{nodal} = d_{proc} + d_{queue} + d_{trans} + d_{prop}$$
$$\text{Traffic Intensity } I = \frac{L \cdot a}{R} \quad (a = \text{average packet arrival rate in packets/sec})$$
$$\text{End-to-End Delay (N links, N-1 routers, ignoring queue/proc)} = N \cdot \left(\frac{L}{R} + \frac{d}{s}\right)$$
$$\text{Throughput} = \min(R_s, R_c, R/\text{connections})$$

---

### [Slide 38] Practice Problem Set on Transmission Delay

#### Problem 1:
- **Question**: A packet of 2,000 bits is transmitted over a 1 Mbps link. Find the transmission delay.
- **Given**:
  - Packet length $L = 2000\text{ bits}$
  - Transmission rate $R = 1\text{ Mbps} = 10^6\text{ bps}$
- **Formula**: $d_{trans} = \frac{L}{R}$
- **Calculation**:
  $$d_{trans} = \frac{2000}{10^6} = 0.002\text{ seconds} = 2\text{ ms}$$
- **Answer**: **$2\text{ ms}$**

#### Problem 2:
- **Question**: A 1,500-byte packet is transmitted over a 100 Mbps Ethernet link. Find the transmission delay.
- **Given**:
  - $L = 1500\text{ bytes} = 1500 \times 8\text{ bits} = 12,000\text{ bits}$
  - $R = 100\text{ Mbps} = 100 \times 10^6\text{ bps} = 10^8\text{ bps}$
- **Calculation**:
  $$d_{trans} = \frac{12000}{10^8} = 0.00012\text{ seconds} = 0.12\text{ ms} = 120\;\mu\text{s}$$
- **Answer**: **$0.12\text{ ms}$ (or $120\;\mu\text{s}$)**

#### Problem 3:
- **Question**: A 5 KB file is sent through a 10 Mbps channel. Find the transmission delay.
- **Given**:
  - $L = 5\text{ KB} = 5 \times 1024 \times 8\text{ bits} = 40,960\text{ bits}$ (or using metric $5 \times 1000 \times 8 = 40,000\text{ bits}$)
  - Using standard networking conversion ($1\text{ KB} = 1024\text{ bytes} = 8192\text{ bits}$):
    $$L = 40,960\text{ bits}$$
  - $R = 10\text{ Mbps} = 10 \times 10^6\text{ bps}$
- **Calculation**:
  $$d_{trans} = \frac{40960}{10^7} = 0.004096\text{ seconds} = 4.096\text{ ms}$$
  *(If using metric $1\text{ KB} = 1000\text{ B}$: $d_{trans} = \frac{40,000}{10^7} = 4\text{ ms}$)*
- **Answer**: **$4.096\text{ ms}$**

#### Problem 4:
- **Question**: A 10 MB file is transmitted over a 50 Mbps connection. Find the transmission delay.
- **Given**:
  - $L = 10\text{ MB} = 10 \times 10^6 \times 8\text{ bits} = 8 \times 10^7\text{ bits}$ (metric) or $10 \times 2^{20} \times 8 = 83,886,080\text{ bits}$
  - $R = 50\text{ Mbps} = 50 \times 10^6\text{ bps}$
- **Calculation**:
  $$d_{trans} = \frac{8 \times 10^7}{50 \times 10^6} = \frac{80}{50} = 1.6\text{ seconds}$$
  *(Binary standard: $\frac{83,886,080}{50,000,000} \approx 1.677\text{ seconds}$)*
- **Answer**: **$1.6\text{ seconds}$**

#### Problem 5:
- **Question**: A transmission delay of 4 ms is observed on a 5 Mbps link. Find the packet size.
- **Given**:
  - $d_{trans} = 4\text{ ms} = 4 \times 10^{-3}\text{ s}$
  - $R = 5\text{ Mbps} = 5 \times 10^6\text{ bps}$
- **Formula**: $L = d_{trans} \times R$
- **Calculation**:
  $$L = (4 \times 10^{-3}\text{ s}) \times (5 \times 10^6\text{ bps}) = 20,000\text{ bits} = 2,500\text{ bytes}$$
- **Answer**: **$20,000\text{ bits}$ ($2.5\text{ KB}$)**

#### Problem 6:
- **Question**: A 12,000-bit packet is transmitted over:
  - Link A = 10 Mbps
  - Link B = 100 Mbps  
  Find the transmission delay for both.
- **Given**: $L = 12,000\text{ bits}$
- **Calculation**:
  - For Link A:
    $$d_{trans, A} = \frac{12000}{10 \times 10^6} = 0.0012\text{ s} = 1.2\text{ ms}$$
  - For Link B:
    $$d_{trans, B} = \frac{12000}{100 \times 10^6} = 0.00012\text{ s} = 0.12\text{ ms} = 120\;\mu\text{s}$$
- **Answer**: **Link A = $1.2\text{ ms}$, Link B = $0.12\text{ ms}$**

---

### [Slide 50] Packet Switching vs. Circuit Switching User Capacity
- **Problem Statement**:
  - Total link capacity $R = 1\text{ Gbps} = 1000\text{ Mbps}$.
  - Each user requires $100\text{ Mbps}$ when active.
  - Each user is active only $10\% = 0.10$ of the time.
  - **Q1**: How many users can use this network under **circuit switching**?
  - **Q2**: How many users can be supported under **packet switching** with near-zero congestion probability?
- **Solution**:
  - **Circuit Switching**:
    $$\text{Number of users } N = \frac{\text{Total Bandwidth}}{\text{Bandwidth per user}} = \frac{1000\text{ Mbps}}{100\text{ Mbps}} = 10\text{ users}$$
    Each user has a reserved circuit. Only **10 users** can be supported simultaneously.
  - **Packet Switching**:
    - With $N = 35$ users, the probability that a given user is active is $p = 0.1$.
    - The probability that more than 10 users are active simultaneously is given by the Binomial distribution:
      $$P(\text{more than 10 active}) = \sum_{k=11}^{35} \binom{35}{k} (0.1)^k (0.9)^{35-k} \approx 0.0004$$
    - Because this probability is less than $0.04\%$, packet switching easily supports **35 users** (more than 3 times circuit switching!) with virtually no noticeable queuing delay.

---

### [Slide 65–66] The Toll Booth / Caravan Analogy Problem
- **Scenario A [Slide 65]**:
  - Caravan of 10 cars (analogous to a 10-bit packet).
  - Toll booth takes 12 seconds to service each car (bit transmission time).
  - Cars propagate at $100\text{ km/h}$ over a $100\text{ km}$ highway.
  - **Q**: How long until the entire caravan is lined up before the second toll booth?
  - **Step 1 (Transmission Delay)**: Time to service/push all 10 cars through the toll booth:
    $$10 \times 12\text{ s} = 120\text{ s} = 2\text{ minutes}$$
  - **Step 2 (Propagation Delay)**: Time for the last car to travel the $100\text{ km}$ distance:
    $$\frac{100\text{ km}}{100\text{ km/h}} = 1\text{ hour} = 60\text{ minutes}$$
  - **Total Time**:
    $$2\text{ minutes} + 60\text{ minutes} = 62\text{ minutes}$$
  - **Answer**: **$62\text{ minutes}$**

- **Scenario B [Slide 66]**:
  - Suppose cars propagate at $1000\text{ km/h}$ and toll booth takes $1\text{ minute}$ per car.
  - **Q**: Will cars arrive at the second booth before all cars are serviced at the first booth?
  - **Step 1**: Time for first car to clear first booth = $1\text{ minute}$.
  - **Step 2**: Propagation time of first car = $\frac{100\text{ km}}{1000\text{ km/h}} = 0.1\text{ hour} = 6\text{ minutes}$.
  - First car arrives at second booth at $t = 1 + 6 = 7\text{ minutes}$.
  - At $t = 7\text{ minutes}$, the first booth has serviced 7 cars. 3 cars are still waiting at the first booth!
  - **Answer**: **Yes, after 7 minutes the first car arrives at the second booth while 3 cars are still at the first booth.**

---

### [Slide 69] Processing Delay and Queueing Delay
- **Problem 1**: A router takes $3\text{ milliseconds}$ to examine every incoming packet before forwarding it. Find the Processing Delay.
  - **Answer**: **$3\text{ ms}$**
- **Problem 2**: A router can transmit one packet every $5\text{ ms}$. When your packet arrives, 4 packets are already waiting in the queue. Calculate the Queueing Delay.
  - **Calculation**:
    $$\text{Queueing Delay} = \text{Number of packets ahead} \times \text{Time per packet} = 4 \times 5\text{ ms} = 20\text{ ms}$$
  - **Answer**: **$20\text{ ms}$**

---

### [Slide 70] Propagation Delay Calculation
- **Given**:
  - Distance $d = 300\text{ km} = 300,000\text{ m} = 3 \times 10^5\text{ m}$
  - Propagation speed $s = 2 \times 10^8\text{ m/s}$
- **Formula**: $d_{prop} = \frac{d}{s}$
- **Calculation**:
  $$d_{prop} = \frac{300,000}{2 \times 10^8} = 1.5 \times 10^{-3}\text{ s} = 1.5\text{ ms}$$
- **Answer**: **$1.5\text{ ms}$**

---

### [Slide 76–78] End-to-End Throughput Calculation
- **Scenario**:
  - Server sending file to client through a router.
  - Server-to-router link capacity: $R_s$
  - Router-to-client link capacity: $R_c$
- **Case 1**: $R_s < R_c$
  - Bottleneck is server link.
  - **Average End-to-End Throughput** = **$R_s\text{ bps}$**
- **Case 2**: $R_s > R_c$
  - Bottleneck is client link.
  - **Average End-to-End Throughput** = **$R_c\text{ bps}$**
- **General Formula for 10 shared flows on core link $R$**:
  $$\text{Throughput} = \min\left(R_s, R_c, \frac{R}{10}\right)$$

---

### [Slide 84–90] Network Topology Formulas
- **Full Mesh Topology**:
  - Number of nodes = $N$
  - Total number of physical duplex links required:
    $$\text{Number of links} = \frac{N(N - 1)}{2}$$
  - Total number of I/O ports required per device:
    $$\text{Ports per device} = N - 1$$
  - Total network ports = $N(N - 1)$
  - *Example*: For $N = 6$ nodes: $\frac{6 \times 5}{2} = 15$ cables, 5 ports per device.
- **Star Topology**:
  - Number of links = $N$
  - Number of ports per node = $1$
  - Number of ports on central hub/switch = $N$
- **Bus Topology**:
  - Number of backbone cables = $1$
  - Number of drop lines = $N$

---

<a name="chapter-2-numericals"></a>
## CHAPTER 2: APPLICATION LAYER NUMERICALS

### [Slide 27–28] HTTP Response Time Calculations

#### Non-Persistent HTTP:
- For each referenced Web object:
  - Step 1: 1 RTT to establish TCP connection (SYN $\to$ SYN-ACK).
  - Step 2: 1 RTT for HTTP request and first bytes of HTTP response.
  - Step 3: File transmission time ($L/R$).
- **Total Non-Persistent HTTP Response Time per object**:
  $$T_{non-persistent} = 2 \cdot RTT + \frac{L}{R}$$
- If a base HTML references $M$ embedded images fetched sequentially:
  $$T_{total} = (2 \cdot RTT + \text{base file trans}) + M \cdot (2 \cdot RTT + \text{image trans})$$

#### Persistent HTTP (HTTP 1.1):
- TCP connection remains open after base HTML is retrieved.
- All $M$ referenced objects are requested back-to-back:
  $$T_{persistent} = 2 \cdot RTT + \text{base file trans} + 1 \cdot RTT + \sum \text{image trans}$$
- Cuts response time nearly in half!

---

### [Slide 47–50] Web Caching (Proxy Server) Numerical Problem

#### Scenario Description:
- Institutional LAN: $1\text{ Gbps}$
- Access link connecting institutional network to public Internet: $R = 1.54\text{ Mbps}$
- Two-way Internet RTT from router to origin server: $RTT = 2.00\text{ seconds}$
- Average Web object size: $L = 100\text{ Kbits} = 100,000\text{ bits}$
- Average request rate from browsers inside institution: $a = 15\text{ requests/sec}$
- Total average data arrival rate requested:
  $$\text{Data rate} = a \times L = 15 \times 100,000\text{ bps} = 1.50\text{ Mbps}$$

#### Calculation Without Cache [Slide 47–48]:
1. **Access Link Utilization**:
   $$\rho_{access} = \frac{\text{Arrival rate}}{\text{Link capacity}} = \frac{1.50\text{ Mbps}}{1.54\text{ Mbps}} \approx 0.97 \quad (97\%)$$
2. **LAN Utilization**:
   $$\rho_{LAN} = \frac{1.50\text{ Mbps}}{1\text{ Gbps}} = 0.0015 \quad (0.15\%)$$
3. **Queueing Delay & Total End-to-End Delay**:
   - Because access link utilization is $0.97$ (approaching 1), queueing delay grows dramatically (several minutes!).
   - End-to-End Delay = $\text{Internet delay (2s)} + \text{Access link delay (minutes)} + \text{LAN delay (microseconds)} \approx \mathbf{\text{several minutes}}$.

#### Calculation With Web Cache [Slide 50]:
- Suppose institutional network installs a local Web cache with a **cache hit rate of 40% ($h = 0.4$)**:
  - $40\%$ of requests are satisfied immediately by local cache (delay $\approx$ tens of milliseconds).
  - $60\%$ ($1 - h = 0.6$) of requests must traverse the access link to origin servers.
1. **New Arrival Rate over Access Link**:
   $$0.60 \times 1.50\text{ Mbps} = 0.90\text{ Mbps}$$
2. **New Access Link Utilization**:
   $$\rho_{access}' = \frac{0.90\text{ Mbps}}{1.54\text{ Mbps}} \approx 0.58 \quad (58\%)$$
   - At $58\%$ utilization, queueing delay drops from minutes to negligible milliseconds!
3. **Average End-to-End Delay**:
   $$\text{Avg Delay} = h \cdot (\text{Cache Delay}) + (1 - h) \cdot (\text{Origin Delay})$$
   $$\text{Avg Delay} = 0.40 \times (0.01\text{ s}) + 0.60 \times (2.01\text{ s}) \approx \mathbf{1.2\text{ seconds}}$$
- **Conclusion**: Installing a cheap local Web cache reduces average user response time from **minutes** down to **$1.2\text{ seconds}$** without purchasing an expensive 100 Mbps leased line.

---

### [Slide 89–92] File Distribution Time: Client-Server vs. P2P

#### Variables:
- $F$: Size of file to distribute (bits)
- $N$: Number of peer clients wanting the file
- $u_s$: Server upload capacity
- $d_{min}$: Minimum client download speed ($d_{min} = \min_i d_i$)
- $u_i$: Upload speed of client peer $i$

#### Mathematical Formulas:
- **Client-Server Distribution Time ($D_{c-s}$)**:
  - Server must upload $N$ copies: time $= \frac{NF}{u_s}$.
  - Slowest client must download one copy: time $= \frac{F}{d_{min}}$.
  $$D_{c-s} \ge \max\left( \frac{NF}{u_s},\; \frac{F}{d_{min}} \right)$$
  *(Increases linearly with $N$ without bound)*

- **P2P Distribution Time ($D_{P2P}$)**:
  - Server must upload at least one copy: $\frac{F}{u_s}$.
  - Slowest client must download one copy: $\frac{F}{d_{min}}$.
  - The aggregate upload capacity of server and all peers is $u_s + \sum_{i=1}^N u_i$.
  $$D_{P2P} \ge \max\left( \frac{F}{u_s},\; \frac{F}{d_{min}},\; \frac{NF}{u_s + \sum_{i=1}^N u_i} \right)$$

#### Worked Example [Slide 92]:
- Given: Client upload rate $u_i = u$, $\frac{F}{u} = 1\text{ hour}$, $u_s = 10u$, $d_{min} \ge u_s$.
- **For $N = 10$ peers**:
  - $D_{c-s} = \max\left(\frac{10F}{10u}, \dots\right) = \frac{F}{u} = 1\text{ hour}$
  - $D_{P2P} = \max\left(\frac{F}{10u}, \frac{10F}{10u + 10u}\right) = \frac{10}{20} \times 1\text{ hr} = 0.5\text{ hour}$
- **For $N = 100$ peers**:
  - $D_{c-s} = \frac{100F}{10u} = 10 \times \frac{F}{u} = \mathbf{10\text{ hours}}$
  - $D_{P2P} = \frac{100F}{10u + 100u} = \frac{100}{110} \times 1\text{ hr} \approx \mathbf{0.91\text{ hour}}$
- **For $N = 1000$ peers**:
  - $D_{c-s} = 100\text{ hours}$
  - $D_{P2P} \approx 1\text{ hour}$ (remains constant and sub-linear!)

---

<a name="chapter-3-numericals"></a>
## CHAPTER 3: TRANSPORT LAYER NUMERICALS

### [Slide 34–37] UDP / Internet Checksum Calculation
- **Algorithm**:
  1. Arrange data into 16-bit integers.
  2. Add them using standard binary addition.
  3. If there is a carry-out from the most significant (17th) bit, wrap around and add 1 to the least significant bit.
  4. Take the 1's complement (invert all bits: $0 \to 1, 1 \to 0$).
- **Worked Slide Example [Slide 36]**:
  - Add two 16-bit integers:
    ```
      1 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0
    + 1 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
    ---------------------------------
    1 1 0 1 1 1 0 1 1 1 0 1 1 1 0 1 1   (17-bit sum with carry-out '1')
    ```
  - **Wraparound (add carry-out to LSB)**:
    ```
      1 0 1 1 1 0 1 1 1 0 1 1 1 1 0 1 1
    +                                 1
    -----------------------------------
      1 0 1 1 1 0 1 1 1 0 1 1 1 1 0 0   (Sum)
    ```
  - **Invert (1's complement)**:
    ```
      0 1 0 0 0 1 0 0 0 1 0 0 0 0 1 1   (Checksum)
    ```
  - **Receiver Verification**:
    - Receiver adds the two numbers plus the checksum:
    $$(1011101110111100)_2 + (0100010001000011)_2 = (1111111111111111)_2$$
    - All 1's confirms no bit error detected.

---

### [Slide 65] Stop-and-Wait Utilization Calculation
- **Problem**:
  - Link speed $R = 1\text{ Gbps} = 10^9\text{ bps}$
  - Propagation delay $d_{prop} = 15\text{ ms}$ (Round Trip Time $RTT = 2 \times 15 = 30\text{ ms}$)
  - Packet size $L = 8,000\text{ bits} = 1\text{ KB}$
- **Calculations**:
  - Packet transmission delay:
    $$d_{trans} = \frac{L}{R} = \frac{8000\text{ bits}}{10^9\text{ bps}} = 8 \times 10^{-6}\text{ s} = 0.008\text{ ms} = 8\;\mu\text{s}$$
  - Sender utilization $U_{sender}$ (fraction of time sender is busy transmitting):
    $$U_{sender} = \frac{L/R}{RTT + L/R} = \frac{0.008\text{ ms}}{30\text{ ms} + 0.008\text{ ms}} = \frac{0.008}{30.008} \approx \mathbf{0.00027 \quad (0.027\%)}$$
  - **Effective Throughput**:
    $$\text{Throughput} = U_{sender} \times R = 0.00027 \times 10^9\text{ bps} = 270\text{ kbps}$$
  - **Conclusion**: A 1 Gbps ultra-high-speed link is utilized at only 270 kbps because Stop-and-Wait protocol requires waiting for an ACK before sending the next packet!

---

### [Slide 69] Pipelining with Window Size $N$
- If sender sends $N$ packets pipelined before waiting for ACK:
  $$U_{sender} = \frac{N \cdot (L/R)}{RTT + (L/R)}$$
- For $N = 3$:
  $$U_{sender} = 3 \times 0.00027 = 0.00081$$
- To achieve $100\%$ link utilization ($U_{sender} = 1$):
  $$N \ge \frac{RTT + L/R}{L/R} = \frac{30.008}{0.008} \approx 3,751\text{ packets}$$

---

### [Slide 77–78] Selective Repeat Window vs. Sequence Number Size
- **Dilemma**: In Selective Repeat, if sequence numbers wrap around too quickly, the receiver cannot distinguish whether an arriving packet is a retransmission of an old packet or a new packet.
- **Formula**:
  $$\text{Window Size } N \le \frac{1}{2} \cdot 2^k = 2^{k-1}$$
  $$\text{Sequence Number Space } \ge 2 \cdot N$$
- *Example from Slide 77*:
  - If sequence numbers are $0, 1, 2, 3$ ($2^2 = 4$ numbers, $k=2$), and window size $N = 3$:
  - $3 > 2^{2-1} = 2$.
  - This causes failure: packets with sequence number 0 arrive, and the receiver cannot tell if it is duplicate packet 0 or new packet 0!
  - Therefore, for window size $N$, sequence number space must be at least $2N$.

---

### [Slide 84–86] TCP Round Trip Time (RTT) and Timeout Estimation
- **Formulas**:
  1. **Estimated RTT (Exponential Weighted Moving Average - EWMA)**:
     $$\text{EstimatedRTT} = (1 - \alpha) \cdot \text{EstimatedRTT} + \alpha \cdot \text{SampleRTT} \quad (\text{typically } \alpha = 0.125)$$
  2. **RTT Deviation (DevRTT)**:
     $$\text{DevRTT} = (1 - \beta) \cdot \text{DevRTT} + \beta \cdot |\text{SampleRTT} - \text{EstimatedRTT}| \quad (\text{typically } \beta = 0.25)$$
  3. **Timeout Interval**:
     $$\text{TimeoutInterval} = \text{EstimatedRTT} + 4 \cdot \text{DevRTT}$$

#### Numerical Example:
- Suppose current $\text{EstimatedRTT} = 100\text{ ms}$, $\text{DevRTT} = 5\text{ ms}$.
- A new packet has $\text{SampleRTT} = 120\text{ ms}$.
- **New EstimatedRTT**:
  $$\text{EstimatedRTT}_{new} = (1 - 0.125) \times 100 + 0.125 \times 120 = 87.5 + 15 = 102.5\text{ ms}$$
- **New DevRTT**:
  $$\text{DevRTT}_{new} = (1 - 0.25) \times 5 + 0.25 \times |120 - 100| = 3.75 + 0.25 \times 20 = 3.75 + 5 = 8.75\text{ ms}$$
- **New TimeoutInterval**:
  $$\text{TimeoutInterval} = 102.5 + 4 \times (8.75) = 102.5 + 35 = \mathbf{137.5\text{ ms}}$$

---

### [Slide 153–154] TCP Throughput & High-Speed Long Fat Pipes
- **Average TCP Throughput Formula (Mathis et al. 1997)**:
  $$\text{Average Throughput} = \frac{1.22 \times MSS}{RTT \times \sqrt{L}}$$
  where $MSS$ is Maximum Segment Size, $RTT$ is round trip time, and $L$ is packet loss probability.
- **Worked Slide Problem [Slide 154]**:
  - $MSS = 1500\text{ bytes} = 12,000\text{ bits}$
  - $RTT = 100\text{ ms} = 0.1\text{ s}$
  - Target throughput = $10\text{ Gbps} = 10^{10}\text{ bps}$
  - **Required Window Size**:
    $$W = \frac{10^{10}\text{ bps} \times 0.1\text{ s}}{12,000\text{ bits}} \approx \mathbf{83,333\text{ segments in flight}}$$
  - **Required Loss Rate $L$**:
    $$10^{10} = \frac{1.22 \times 12000}{0.1 \times \sqrt{L}} \implies \sqrt{L} = \frac{14,640}{10^9} = 1.464 \times 10^{-5}$$
    $$L = (1.464 \times 10^{-5})^2 \approx \mathbf{2 \times 10^{-10}}$$
  - **Takeaway**: To sustain 10 Gbps TCP throughput over a 100ms link, packet loss probability cannot exceed 1 packet in every 5 billion packets!

---

<a name="chapter-4-numericals"></a>
## CHAPTER 4: NETWORK LAYER (DATA PLANE) NUMERICALS

### [Slide 33] Router Buffer Sizing Rule of Thumb
- **Traditional RFC 3439 Rule of Thumb**:
  $$\text{Buffer Size } B = RTT \times C$$
  where $RTT \approx 250\text{ ms}$, $C = \text{link capacity}$.
  - *Example*: For a $C = 10\text{ Gbps}$ link:
    $$B = 0.25\text{ s} \times 10 \times 10^9\text{ bps} = \mathbf{2.5\text{ Gbits buffer}}$$
- **Modern Rule for $N$ Independent TCP Flows**:
  $$B = \frac{RTT \times C}{\sqrt{N}}$$
  - *Example*: If $N = 100$ flows on a 10 Gbps link:
    $$B = \frac{2.5\text{ Gbits}}{\sqrt{100}} = \frac{2.5\text{ Gbits}}{10} = \mathbf{250\text{ Mbits}}$$

---

### [Slide 54–62] IP Addressing, Subnetting & CIDR Calculations

#### Classful Addressing Summary Table [Slides 55–59]:
| Class | Leading Bits | Network ID Bits | Host ID Bits | Default Subnet Mask | Total Networks | Usable Hosts / Network ($2^h - 2$) | IP Range |
|---|---|---|---|---|---|---|---|
| **A** | `0` | 8 (7 usable) | 24 | `255.0.0.0` (/8) | $2^7 = 128$ | $2^{24} - 2 = 16,777,214$ | `0.0.0.0` - `127.255.255.255` |
| **B** | `10` | 16 (14 usable)| 16 | `255.255.0.0` (/16) | $2^{14} = 16,384$ | $2^{16} - 2 = 65,534$ | `128.0.0.0` - `191.255.255.255` |
| **C** | `110` | 24 (21 usable)| 8 | `255.255.255.0` (/24) | $2^{21} = 2,097,152$ | $2^8 - 2 = 254$ | `192.0.0.0` - `223.255.255.255` |
| **D** | `1110` | Multicast (28 bits group) | None | N/A | N/A | N/A | `224.0.0.0` - `239.255.255.255` |
| **E** | `1111` | Reserved / Experimental | None | N/A | N/A | N/A | `240.0.0.0` - `255.255.255.255` |

#### Special IP Addresses [Slide 60]:
- `127.0.0.0/8`: Loopback / Localhost (tests TCP/IP stack locally).
- `169.254.0.0/16`: APIPA (Automatic Private IP Addressing) / Link-Local when DHCP fails.
- `0.0.0.0/8`: Represents current network during DHCP initialization.
- Private IP Address Blocks (RFC 1918):
  - `10.0.0.0/8` ($10.0.0.0$ - $10.255.255.255$)
  - `172.16.0.0/12` ($172.16.0.0$ - $172.31.255.255$)
  - `192.168.0.0/16` ($192.168.0.0$ - $192.168.255.255$)

#### Subnetting & CIDR Address Allocation [Slide 71]:
- ISP owns block: `200.23.16.0/20` (`11001000 00010111 00010000 00000000`)
- Total addresses in /20 = $2^{32-20} = 2^{12} = 4096$ addresses.
- ISP divides its block into 8 equal subnets for 8 organizations:
  - Each organization gets $/23$ ($20 + 3$ bits borrowed $\implies 2^3 = 8$ subnets):
  - Addresses per org: $2^{32-23} = 2^9 = 512$ addresses ($510$ usable hosts).
  - **Subnet Allocations**:
    - Org 0: `200.23.16.0/23` (`200.23.16.0` - `200.23.17.255`)
    - Org 1: `200.23.18.0/23` (`200.23.18.0` - `200.23.19.255`)
    - Org 2: `200.23.20.0/23` (`200.23.20.0` - `200.23.21.255`)
    - Org 3: `200.23.22.0/23` (`200.23.22.0` - `200.23.23.255`)
    - Org 4: `200.23.24.0/23` (`200.23.24.0` - `200.23.25.255`)
    - Org 5: `200.23.26.0/23` (`200.23.26.0` - `200.23.27.255`)
    - Org 6: `200.23.28.0/23` (`200.23.28.0` - `200.23.29.255`)
    - Org 7: `200.23.30.0/23` (`200.23.30.0` - `200.23.31.255`)

---

### [Slide 84–88] NAT Translation Table Numerical Examples
- **Example 1: Basic Outbound HTTP [Slide 85]**:
  - Internal Host `10.10.5.20`, Source Port `49200` connects to Web Server `198.51.100.20:80`.
  - Public NAT IP = `203.0.113.9`.
  - NAT allocates WAN port `30001`.
  - **NAT Table Entry**: `203.0.113.9:30001` $\longleftrightarrow$ `10.10.5.20:49200`.
  - Outgoing packet rewritten to: Source `203.0.113.9:30001`, Dest `198.51.100.20:80`.
  - Return packet: Source `198.51.100.20:80`, Dest `203.0.113.9:30001`.
  - NAT table lookup rewrites Dest to `10.10.5.20:49200`.

- **Example 2: Same Internal Port Collision [Slide 87]**:
  - Laptop `192.168.0.5` and Phone `192.168.0.6` both pick ephemeral port `50000`.
  - NAT Router IP: `203.0.113.77`.
  - Laptop sends first: NAT maps `192.168.0.5:50000` $\longleftrightarrow$ `203.0.113.77:50000`.
  - Phone sends next: Public port 50000 is taken! NAT assigns new port `50001`.
  - Phone mapped to: `192.168.0.6:50000` $\longleftrightarrow$ `203.0.113.77:50001`.
  - Distinct public port numbers allow a single public IP to multiplex thousands of internal hosts!

---

### [Slide 122] IPv4 Datagram Fragmentation and Reassembly
- **Problem**:
  - Original Datagram: Total length = $4000\text{ bytes}$ (includes $20\text{ bytes}$ IP header and $3980\text{ bytes}$ data payload).
  - Outgoing Link MTU (Maximum Transmission Unit) = $1500\text{ bytes}$.
  - Fragment data size must be a multiple of 8 bytes.
  - $\max \text{ payload per fragment} = 1500 - 20 = 1480\text{ bytes}$ (since $1480 / 8 = 185$, valid multiple of 8).
- **Fragment Calculations**:
  - **Fragment 1**:
    - Total Length = $1480 + 20 = 1500\text{ bytes}$
    - Identification = $x$
    - Flags: More Fragments (MF) = `1` (more fragments follow)
    - Fragmentation Offset = $\frac{0}{8} = \mathbf{0}$
    - Data carried: bytes $0$ to $1479$ ($1480\text{ bytes}$)
  - **Fragment 2**:
    - Total Length = $1480 + 20 = 1500\text{ bytes}$
    - Identification = $x$
    - Flags: More Fragments (MF) = `1`
    - Fragmentation Offset = $\frac{1480}{8} = \mathbf{185}$
    - Data carried: bytes $1480$ to $2959$ ($1480\text{ bytes}$)
  - **Fragment 3**:
    - Remaining payload = $3980 - 1480 - 1480 = 1020\text{ bytes}$.
    - Total Length = $1020 + 20 = 1040\text{ bytes}$
    - Identification = $x$
    - Flags: More Fragments (MF) = `0` (last fragment!)
    - Fragmentation Offset = $\frac{1480 + 1480}{8} = \frac{2960}{8} = \mathbf{370}$
    - Data carried: bytes $2960$ to $3979$ ($1020\text{ bytes}$)

---

<a name="chapter-5-numericals"></a>
## CHAPTER 5: NETWORK LAYER (CONTROL PLANE) NUMERICALS

### [Slide 14–27] Dijkstra's Link-State Routing Algorithm
- **Network Topology**:
  - Nodes: $u, v, w, x, y, z$
  - Source Node: $u$
  - Direct Link Costs:
    - $c(u,v) = 2$, $c(u,x) = 1$, $c(u,w) = 5$
    - $c(x,v) = 2$, $c(x,w) = 3$, $c(x,y) = 1$
    - $c(v,w) = 3$, $c(w,y) = 1$, $c(w,z) = 5$
    - $c(y,z) = 2$
- **Algorithm Execution Table**:

| Step | Set $N'$ | $D(v), p(v)$ | $D(w), p(w)$ | $D(x), p(x)$ | $D(y), p(y)$ | $D(z), p(z)$ | Chosen Node (min) |
|---|---|---|---|---|---|---|---|
| **0** | $\{u\}$ | $2, u$ | $5, u$ | **$1, u$** | $\infty$ | $\infty$ | Node $x$ |
| **1** | $\{u, x\}$ | $2, u$ | $\min(5, 1+3) = 4, x$ | --- | **$\min(\infty, 1+1) = 2, x$** | $\infty$ | Node $y$ (tie with $v$) |
| **2** | $\{u, x, y\}$ | **$2, u$** | $\min(4, 2+1) = 3, y$ | --- | --- | $\min(\infty, 2+2) = 4, y$ | Node $v$ |
| **3** | $\{u, x, y, v\}$ | --- | **$3, y$** | --- | --- | $4, y$ | Node $w$ |
| **4** | $\{u, x, y, v, w\}$| --- | --- | --- | --- | **$4, y$** | Node $z$ |
| **5** | $\{u, x, y, v, w, z\}$| Full tree constructed | | | | | Finished |

- **Shortest Paths from $u$**:
  - To $x$: Path $u-x$, cost = **1**
  - To $v$: Path $u-v$, cost = **2**
  - To $y$: Path $u-x-y$, cost = **2**
  - To $w$: Path $u-x-y-w$, cost = **3**
  - To $z$: Path $u-x-y-z$, cost = **4**
- **Forwarding Table at router $u$**:
  - Destination $x$: outgoing link $(u,x)$
  - Destination $v$: outgoing link $(u,v)$
  - Destination $y$: outgoing link $(u,x)$
  - Destination $w$: outgoing link $(u,x)$
  - Destination $z$: outgoing link $(u,x)$

---

### [Slide 34–50, 54–56] Bellman-Ford Distance Vector Algorithm
- **Equation**:
  $$d_x(y) = \min_v \{ c(x,v) + d_v(y) \}$$
- **Initial Matrix ($t=0$) [Slide 54]**:

| Node | $u$ | $v$ | $x$ | $w$ | $y$ | $z$ |
|---|---|---|---|---|---|---|
| **$u$** | 0 | 2 | 1 | 5 | $\infty$ | $\infty$ |
| **$v$** | 2 | 0 | 2 | 3 | $\infty$ | $\infty$ |
| **$x$** | 1 | 2 | 0 | 3 | 1 | $\infty$ |
| **$w$** | 5 | 3 | 3 | 0 | 1 | 5 |
| **$y$** | $\infty$ | $\infty$ | 1 | 1 | 0 | 2 |
| **$z$** | $\infty$ | $\infty$ | $\infty$ | 5 | 2 | 0 |

- **Iteration 1 & Calculation for $u$**:
  - $D_u(y) = \min\{c(u,v)+D_v(y),\; c(u,x)+D_x(y),\; c(u,w)+D_w(y)\} = \min\{2+\infty,\; 1+1,\; 5+1\} = 2$ (via $x$).
  - $D_u(z) = \min\{c(u,v)+D_v(z),\; c(u,x)+D_x(z),\; c(u,w)+D_w(z)\} = \min\{2+\infty,\; 1+\infty,\; 5+5\} = 10$ (via $w$).
- **Iteration 2 (Convergence)**:
  - After neighbor updates propagate:
  - $D_u(z) = \min\{c(u,x) + D_x(z)\} = 1 + (D_x(y) + D_y(z)) = 1 + (1 + 2) = 4$ (via $x$).
  - Converged minimum costs: $u \to x = 1, u \to v = 2, u \to y = 2, u \to w = 3, u \to z = 4$.

---

### [Slide 52–53] Count-to-Infinity Problem
- **Scenario**:
  - Link between $x$ and $y$ increases from cost 4 to cost 60.
  - $z$ had route to $x$ via $y$ of cost 5.
  - $y$ detects link cost change to 60. It sees $z$ advertised path to $x$ of cost 5.
  - $y$ incorrectly thinks: "I can reach $x$ through $z$ with cost $c(y,z) + D_z(x) = 1 + 5 = 6$".
  - $y$ tells $z$: $D_y(x) = 6$.
  - $z$ updates: $D_z(x) = c(z,y) + D_y(x) = 1 + 6 = 7$.
  - $y$ updates: $D_y(x) = 1 + 7 = 8$.
  - Values loop and increment $6 \to 7 \to 8 \to 9 \dots$ up to 60 before stabilizing!
- **Solution**: **Poisoned Reverse** (if $z$ routes through $y$ to $x$, $z$ tells $y$ that $D_z(x) = \infty$).

---

<a name="chapter-6-numericals"></a>
## CHAPTER 6: LINK LAYER & LANS NUMERICALS

### [Slide 13] Parity Checking (1D and 2D)
- **Even Parity**: Parity bit is set such that the total number of 1s in the data plus parity bit is **even**.
  - Example: Data = `01110001` (four 1s). Parity bit = `0`.
  - Example: Data = `01110011` (five 1s). Parity bit = `1`.
- **Two-Dimensional Parity**:
  - Calculates row parity and column parity.
  - Can **detect and correct** any single bit error (intersection of invalid row and column identifies the flipped bit).

---

### [Slide 15–18] CRC (Cyclic Redundancy Check) Division
- **Given from Slide 17**:
  - Data bits $D = 101110$ ($6\text{ bits}$)
  - Generator polynomial $G = 1001$ ($r+1 = 4\text{ bits} \implies r = 3\text{ CRC bits}$)
- **Step 1: Append $r = 3$ zeros to $D$**:
  $$D \cdot 2^r = 101110000$$
- **Step 2: Modulo-2 Division (XOR division)**:
  ```
               1 0 1 0 1 1  (Quotient)
         ------------------
  1 0 0 1 ) 1 0 1 1 1 0 0 0 0
            1 0 0 1
            -------
              0 1 0 1
              0 0 0 0
              -------
                1 0 1 0
                1 0 0 1
                -------
                  0 1 1 0
                  0 0 0 0
                  -------
                    1 1 0 0
                    1 0 0 1
                    -------
                      1 0 1 0
                      1 0 0 1
                      -------
                        0 1 1  <--- Remainder R (3 bits)
  ```
- **Step 3: Remainder $R = 011$**:
- **Transmitted Frame**:
  $$\langle D, R \rangle = 101110011$$
- **Verification at Receiver**:
  - Receiver divides $101110011$ by $1001$. The remainder is $000$, confirming no error!

---

### [Slide 29–30, 113] Multiple Access Protocol Efficiencies

#### Slotted ALOHA Efficiency [Slide 29]:
- $N$ active nodes, each transmits in a slot with probability $p$.
- Probability that a given node succeeds: $P(\text{success}) = p (1 - p)^{N-1}$
- Total efficiency (probability that any node succeeds):
  $$E(p) = N \cdot p (1 - p)^{N-1}$$
- Maximum efficiency occurs at $p^* = 1/N$. Taking the limit as $N \to \infty$:
  $$\lim_{N \to \infty} N \cdot \frac{1}{N} \left(1 - \frac{1}{N}\right)^{N-1} = \frac{1}{e} \approx \mathbf{0.37 \quad (37\%)}$$

#### Pure (Unslotted) ALOHA Efficiency [Slide 30, 113]:
- A frame sent at $t_0$ collides with frames sent in $[t_0 - 1, t_0 + 1]$ (vulnerable window is $2 \times$ frame time).
- Probability of no collision = $(1 - p)^{2(N-1)}$.
- Efficiency:
  $$E_{pure} = N \cdot p (1 - p)^{2(N-1)} \xrightarrow{N \to \infty} \frac{1}{2e} \approx \mathbf{0.18 \quad (18\%)}$$

#### CSMA/CD Efficiency [Slide 35]:
$$\text{Efficiency} = \frac{1}{1 + 5 \cdot \frac{t_{prop}}{t_{trans}}}$$
- As propagation delay $t_{prop} \to 0$, Efficiency $\to 100\%$.
- As frame transmission time $t_{trans} \to \infty$, Efficiency $\to 100\%$.

---

### [Slide 71–72] Self-Learning Switch Table Construction
- **Problem**:
  - Host C sends a frame to Host I on switch S1.
  - Switch S1 has empty table.
- **Action**:
  1. S1 learns: `MAC(C)` is on `Port 1`, with TTL.
  2. S1 does not know where `MAC(I)` is $\implies$ S1 **floods** the frame out all ports except Port 1.
  3. Frame reaches S4, S2, S3 which also flood.
  4. When Host I receives the frame, Host I replies with a frame destined to C.
  5. S3 learns `MAC(I)` is on `Port 4`. S3 already knows C is reachable via S4 $\implies$ S3 **forwards selectively** without flooding!

---

<a name="chapter-7-numericals"></a>
## CHAPTER 7: WIRELESS & MOBILE NETWORKS NUMERICALS

### [Slide 20] Signal-to-Noise Ratio (SNR) in Decibels
$$\text{SNR}_{dB} = 10 \log_{10}\left(\frac{S}{N}\right)$$
- If Signal Power = $1000\text{ mW}$ and Noise Power = $1\text{ mW}$:
  $$\frac{S}{N} = 1000 \implies \text{SNR}_{dB} = 10 \log_{10}(1000) = 10 \times 3 = \mathbf{30\text{ dB}}$$
- If $\text{SNR}_{dB} = 0\text{ dB}$: Signal power equals Noise power ($S = N$).

---

### [Slide 21] Shannon Channel Capacity
$$C = B \cdot \log_2\left(1 + \frac{S}{N}\right) \quad \text{(in bits/sec)}$$
where $B$ is channel bandwidth in Hertz, and $S/N$ is linear signal-to-noise ratio.
- **Numerical Example**:
  - Bandwidth $B = 20\text{ MHz} = 20 \times 10^6\text{ Hz}$
  - $\text{SNR} = 30\text{ dB} \implies \frac{S}{N} = 10^{30/10} = 1000$.
  - Maximum Capacity:
    $$C = 20 \times 10^6 \times \log_2(1 + 1000) \approx 20 \times 10^6 \times 9.967 \approx \mathbf{199.3\text{ Mbps}}$$

---

### [Slide 44] Digital Modulation Bits per Symbol
- **BPSK** (Binary Phase Shift Keying): 2 phases $\implies 1\text{ bit/symbol}$
- **QPSK** (Quadrature Phase Shift Keying): 4 phases $\implies \log_2(4) = 2\text{ bits/symbol}$
- **16-QAM**: 16 constellation points $\implies \log_2(16) = 4\text{ bits/symbol}$
- **64-QAM**: 64 constellation points $\implies \log_2(64) = 6\text{ bits/symbol}$
- **256-QAM**: 256 constellation points $\implies \log_2(256) = 8\text{ bits/symbol}$
- *Note*: Higher modulation provides higher data rate but requires higher SNR to maintain low Bit Error Rate (BER).

---

### [Slide 79] 4G/5G Resource Block (RB) Sizing
- A Resource Block consists of:
  - 12 subcarriers (each $15\text{ kHz}$ wide $\implies 180\text{ kHz}$ total bandwidth)
  - 1 time slot of $0.5\text{ ms}$ containing 7 OFDM symbols.
- **Total Resource Elements per RB**:
  $$12 \text{ subcarriers} \times 7 \text{ OFDM symbols} = \mathbf{84\text{ Resource Elements (RE)}}$$

---

<a name="chapter-8-numericals"></a>
## CHAPTER 8: NETWORK SECURITY NUMERICALS

### [Slide 19] Modular Arithmetic Rules
- $(a + b) \bmod n = [(a \bmod n) + (b \bmod n)] \bmod n$
- $(a \times b) \bmod n = [(a \bmod n) \times (b \bmod n)] \bmod n$
- $a^d \bmod n = (a \bmod n)^d \bmod n$
- *Example from Slide 19*:
  - $x = 14, n = 10, d = 2$:
  - $(14 \bmod 10)^2 \bmod 10 = 4^2 \bmod 10 = 16 \bmod 10 = \mathbf{6}$
  - $14^2 \bmod 10 = 196 \bmod 10 = \mathbf{6}$

---

### [Slide 21–24] RSA Public-Key Cryptography Complete Example

#### Key Generation Algorithm [Slide 21]:
1. Choose two prime numbers: $p = 5$, $q = 7$.
2. Compute modulus $n$:
   $$n = p \times q = 5 \times 7 = 35$$
3. Compute Euler's totient $z$:
   $$z = (p - 1)(q - 1) = (5 - 1)(7 - 1) = 4 \times 6 = 24$$
4. Choose public exponent $e$ such that $\gcd(e, z) = 1$ and $e < n$:
   - Let $e = 5$ (since factors of 24 are 2 and 3; 5 is relatively prime to 24).
5. Compute private key $d$ such that $(e \cdot d) \bmod z = 1$:
   $$(5 \times d) \bmod 24 = 1$$
   - Testing values:
     - $d = 5 \implies 25 \bmod 24 = 1$. (Here $d=5$, but to keep $e \ne d$ we can use $d = 29$, since $5 \times 29 = 145 = 6 \times 24 + 1$).
   - Slide chooses **$d = 29$**.
6. **Keys Generated**:
   - **Public Key**: $(n, e) = (35, 5)$
   - **Private Key**: $(n, d) = (35, 29)$

#### Encryption [Slide 23]:
- Message $m = 12$ ($m < n$)
- Ciphertext $c = m^e \bmod n$:
  $$c = 12^5 \bmod 35$$
  - $12^2 = 144 \equiv 4 \pmod{35}$
  - $12^4 = 4^2 = 16 \pmod{35}$
  - $12^5 = 12 \times 16 = 192 \pmod{35}$
  - $192 = 5 \times 35 + 17 \implies c = \mathbf{17}$
- **Ciphertext = 17**

#### Decryption [Slide 23]:
- Plaintext $m = c^d \bmod n$:
  $$m = 17^{29} \bmod 35$$
  - Using modular exponentiation:
  - $17^1 \equiv 17 \pmod{35}$
  - $17^2 = 289 = 8 \times 35 + 9 \equiv 9 \pmod{35}$
  - $17^4 = 9^2 = 81 \equiv 11 \pmod{35}$
  - $17^8 = 11^2 = 121 \equiv 16 \pmod{35}$
  - $17^{16} = 16^2 = 256 \equiv 11 \pmod{35}$
  - $29 = 16 + 8 + 4 + 1$:
  - $17^{29} = 17^{16} \times 17^8 \times 17^4 \times 17^1 \equiv 11 \times 16 \times 11 \times 17 \pmod{35}$
  - $11 \times 16 = 176 \equiv 1 \pmod{35}$
  - $1 \times 11 \times 17 = 187 = 5 \times 35 + 12 \equiv \mathbf{12} \pmod{35}$
- **Recovered Plaintext = 12** (Decryption matches original message!).

---

### [Slide 46] Checksum Hash Collision Demonstration
- Demonstrates why the Internet Checksum cannot be used as a cryptographic hash function.
- **Message 1**: `I O U 1 0 0 . 9 9 B O B`
  - ASCII Hex: `49 4F 55 31 30 30 2E 39 39 42 D2 42`
  - 16-bit Checksum = `B2 C1 D2 AC`
- **Message 2 (Attacker swap)**: `I O U 9 0 0 . 1 9 B O B`
  - Swapped the '1' and '9' in the bill amount!
  - 16-bit Checksum = `B2 C1 D2 AC`
- **Result**: Different messages produce the exact same checksum! A cryptographic hash must be collision-resistant (e.g., SHA-256).

---

# PART 2: THEORETICAL QUESTIONS, "WHY", AND COMPARISONS

<a name="chapter-1-theory"></a>
## CHAPTER 1: FOUNDATIONS & LAYERED ARCHITECTURES

### Q1: What is the difference between Packet Switching and Circuit Switching? [Slides 45–51]
| Feature | Circuit Switching | Packet Switching |
|---|---|---|
| **Resource Allocation** | Dedicated end-to-end circuit reserved for call duration | On-demand resource sharing (statistical multiplexing) |
| **Performance** | Guaranteed circuit-like performance (constant rate) | Variable performance; subject to queuing delay & packet loss |
| **Idle Resource Usage** | Dedicated link capacity wasted when silent | Idle capacity used immediately by other packets |
| **Setup Phase** | Call setup required before data transfer | No call setup needed; packets sent immediately |
| **Best Suited For** | Continuous, real-time voice traffic | Bursty data traffic (Web, file transfer, apps) |

### Q2: Why is Packet Switching well-suited for "bursty" data traffic? [Slide 51]
- Users browsing the web or sending messages only transmit for brief moments and spend most time reading (idle).
- Circuit switching would reserve 100% of bandwidth for an idle user, wasting 90%+ of channel capacity.
- Packet switching allows dozens of users to dynamically share the same link, achieving over 3x the user capacity with negligible probability of buffer overflow.

### Q3: What are the Four Sources of Packet Delay? [Slides 63–64]
1. **Nodal Processing Delay ($d_{proc}$)**: Time required by the router to examine the packet header, check for bit-level errors, and determine the outgoing interface (typically $< \mu\text{s}$).
2. **Queueing Delay ($d_{queue}$)**: Time the packet spends waiting in router memory buffers until the output link becomes free (depends heavily on network congestion).
3. **Transmission Delay ($d_{trans} = L/R$)**: Time required to push all $L$ bits of the packet onto the physical communication link at transmission speed $R$.
4. **Propagation Delay ($d_{prop} = d/s$)**: Time required for a single bit to travel across the physical medium over distance $d$ at physical speed $s$ ($pprox 2 \times 10^8\text{ m/s}$).

### Q4: Compare Network Topologies: Star, Bus, Ring, Mesh, Tree, and Hybrid [Slides 79–96]
| Topology | Cable Requirement | Fault Tolerance / Reliability | Scalability & Cost | Disadvantages |
|---|---|---|---|---|
| **Mesh** | $N(N-1)/2$ cables | Extremely high; alternate paths exist | Very expensive; difficult installation | High cable cost; massive wiring |
| **Star** | $N$ cables | High; single cable failure affects only 1 host | Highly scalable; easy to add nodes | Single point of failure at central hub |
| **Bus** | 1 backbone + $N$ drops | Low; backbone break halts whole network | Cheap; limited cable | Heavy collisions; hard fault troubleshooting |
| **Ring** | $N$ cables (closed ring) | Low (unidirectional ring breaks if 1 node dies)| Moderate | Token delay; single break halts communication |
| **Tree** | Hierarchical star | Good; localized segment failure | Easy to expand in tiers | Failure of root/hub disrupts branch |
| **Hybrid** | Combines 2+ topologies | Flexible | High cost; complex management | High architecture complexity |

### Q5: What is the difference between Physical Topology and Logical Topology? [Slide 79]
- **Physical Topology**: Refers to the physical layout and arrangement of cables, computers, and networking hardware.
- **Logical Topology**: Refers to the logical path and protocol rules data follows to travel across the network (e.g., Ethernet physically wired as a star switch operates logically as a broadcast bus or point-to-point).

### Q6: What are the differences between PAN, LAN, CAN, MAN, and WAN? [Slides 99–103]
| Network Type | Full Form | Typical Geographic Range | Example Technologies |
|---|---|---|---|
| **PAN** | Personal Area Network | $1 - 10\text{ meters}$ | Bluetooth, Zigbee, NFC |
| **LAN** | Local Area Network | Room, building, campus ($100\text{m} - 1\text{km}$)| Ethernet, Wi-Fi |
| **CAN** | Campus Area Network | University or corporate complex ($1 - 5\text{km}$)| Interconnected LANs via fiber |
| **MAN** | Metropolitan Area Network| Entire city ($10 - 50\text{km}$) | Cable TV broadband, Metro Ethernet |
| **WAN** | Wide Area Network | Country, continent, global ($> 100\text{km}$) | The Internet, undersea fiber links |

### Q7: Why Layering? What are the Advantages and Disadvantages? [Slide 142]
- **Why Layering**: Divides a complex communication system into modular, well-defined functional layers with explicit abstractions and interfaces.
- **Advantages**:
  1. Modularity makes system design, maintenance, and debugging simpler.
  2. Modifying an implementation at one layer does not affect other layers (e.g., switching from Wi-Fi to Ethernet does not break HTTP).
- **Disadvantages**:
  1. Potential performance overhead due to duplicate functions (e.g., error checking at link, transport, and application layers).
  2. Information hiding prevents lower layers from utilizing upper layer information (violating cross-layer optimization).

### Q8: Compare the OSI 7-Layer Model and the TCP/IP 5-Layer Stack [Slides 120–138]
- **The OSI 7 Layers**: Physical $\to$ Data Link $\to$ Network $\to$ Transport $\to$ Session $\to$ Presentation $\to$ Application.
- **TCP/IP 5-Layer Stack**: Physical $\to$ Link $\to$ Network (IP) $\to$ Transport (TCP/UDP) $\to$ Application.
- **Two Layers in OSI Not Present in TCP/IP**:
  1. **Presentation Layer (Layer 6)**: Responsible for data formatting, translation (ASCII/EBCDIC), compression, and encryption. In TCP/IP, this is handled by the application developer.
  2. **Session Layer (Layer 5)**: Responsible for session checkpointing, recovery, and synchronization. In TCP/IP, handled within the application layer.

---

<a name="chapter-2-theory"></a>
## CHAPTER 2: APPLICATION LAYER

### Q9: Why does UDP exist when TCP provides reliability? [Slide 17]
- TCP adds substantial overhead: 3-way connection handshake, retransmission timers, ordered buffering, flow control, and congestion throttling.
- UDP provides raw, low-overhead datagram transmission:
  1. **No connection setup delay**: Crucial for DNS queries and real-time gaming.
  2. **No rate throttling**: Streaming audio/video and VoIP can send at constant bitrates without being halved by TCP congestion control.
  3. **Small packet header**: UDP header is 8 bytes vs. TCP's 20 bytes.

### Q10: Difference between Non-Persistent HTTP and Persistent HTTP [Slides 24–28]
- **Non-Persistent (HTTP 1.0)**:
  - At most one object transmitted per TCP connection.
  - Requires $2 \cdot RTT$ per object download.
  - Causes significant OS socket overhead and slow start latency.
- **Persistent (HTTP 1.1)**:
  - TCP connection remains open across multiple object requests.
  - Successive requests sent back-to-back over the open connection.
  - Reduces download time to as little as $1 \cdot RTT$ for all referenced objects.

### Q11: How do Cookies maintain state in a stateless HTTP protocol? [Slides 35–36]
- HTTP is inherently stateless (each GET request is independent).
- **Four Components of Cookies**:
  1. `Set-cookie:` response header sent by Web server.
  2. `Cookie:` request header included automatically by client browser on subsequent requests.
  3. Client cookie file stored on host machine.
  4. Backend database associating cookie ID with shopping carts, login sessions, and preferences.

### Q12: Why not centralize DNS? [Slide 72]
A single centralized DNS server would fail because:
1. **Single point of failure**: If it crashes, the entire global Internet stops working.
2. **Traffic volume**: A single site cannot handle trillions of queries per day (Comcast alone handles 600B/day).
3. **Distant centralized database**: Delays for users across continents would be unacceptable.
4. **Maintenance**: Maintaining records for billions of hosts in one place is unmanageable.
- *Solution*: Distributed, hierarchical database (Root servers $\to$ TLD servers $\to$ Authoritative servers).

### Q13: Difference between Iterated Query and Recursive Query in DNS [Slides 79–80]
- **Iterated Query**: Contacted server replies with the referral address of the next DNS server to ask ("I don't know, but ask this server"). The local DNS server performs all sequential queries.
- **Recursive Query**: Puts the burden of resolution on the contacted server (each contacted server queries the next on behalf of the client). Heavy load on upper root/TLD servers, hence rarely used at the top levels.

### Q14: Explain DNS Resource Record (RR) Types [Slide 82]
Format: `(Name, Value, Type, TTL)`
- `Type=A`: Name is hostname, Value is 32-bit IPv4 address (e.g., `ibm.com` $\to$ `192.0.2.1`).
- `Type=NS`: Name is domain, Value is authoritative name server hostname (e.g., `foo.com` $\to$ `dns.foo.com`).
- `Type=CNAME`: Name is alias, Value is canonical (true) hostname.
- `Type=MX`: Name is domain, Value is mail server hostname.

---

<a name="chapter-3-theory"></a>
## CHAPTER 3: TRANSPORT LAYER

### Q15: Compare Go-Back-N (GBN) and Selective Repeat (SR) [Slides 70–78]
| Feature | Go-Back-N (GBN) | Selective Repeat (SR) |
|---|---|---|
| **Acknowledgements** | Cumulative ACKs (ACKs all packets up to $n$) | Individual ACKs for every correctly received packet |
| **Out-of-order Packets** | Discarded by receiver (not buffered) | Buffered by receiver until gaps are filled |
| **Sender Timers** | Single timer for oldest unACKed packet | Individual logical timer for every unACKed packet |
| **Retransmissions** | Retransmits packet $n$ and **all** subsequent packets | Retransmits **only** the specific packet that timed out |
| **Window Size Limit** | $W \le 2^k - 1$ | $W \le 2^{k-1}$ |

### Q16: Difference between Flow Control and Congestion Control [Slides 93, 109]
- **Flow Control**: Point-to-point mechanism to prevent the sender from overflowing the **receiver's application buffer**. Controlled via the `Receive Window (rwnd)` header field.
- **Congestion Control**: Network-wide mechanism to prevent sending hosts from overloading the **network routers and switches**. Controlled via the `Congestion Window (cwnd)` based on perceived packet loss/delay.

### Q17: Explain TCP 3-Way Handshake and Connection Teardown [Slides 99–104]
- **Connection Establishment (3-Way Handshake)**:
  1. Client sends `SYN` (seq = $x$, ACK = 0).
  2. Server responds with `SYN-ACK` (seq = $y$, ACK = $x+1$).
  3. Client replies with `ACK` (seq = $x+1$, ACK = $y+1$, can carry data).
- **Connection Teardown**:
  1. Initiator sends `FIN`.
  2. Receiver responds with `ACK`.
  3. Receiver sends its own `FIN`.
  4. Initiator replies with `ACK` and waits $2 \times \text{MSL}$ (Maximum Segment Lifetime) in `TIME_WAIT` state before closing.

### Q18: Explain TCP AIMD (Additive Increase Multiplicative Decrease) [Slides 123–126]
- **Additive Increase**: Sender increases `cwnd` by 1 MSS each RTT in the absence of loss (probing for extra bandwidth).
- **Multiplicative Decrease**: When packet loss occurs:
  - On triple duplicate ACK (TCP Reno): `cwnd` is cut in half.
  - On timeout: `cwnd` is reset to 1 MSS, entering Slow Start.

---

<a name="chapter-4-theory"></a>
## CHAPTER 4: NETWORK LAYER (DATA PLANE)

### Q19: What is the difference between Routing and Forwarding? [Slide 42]
- **Forwarding (Data Plane)**: Local per-router action that moves arriving packets from an input port to the appropriate output port based on the forwarding table (timescale: nanoseconds, hardware implemented).
- **Routing (Control Plane)**: Network-wide logic that determines the end-to-end paths taken by packets from source to destination across multiple routers (timescale: seconds, software implemented).

### Q20: Explain Head-of-the-Line (HOL) Blocking in Router Input Ports [Slide 30]
- Occurs when multiple packets are queued at an input port buffer.
- If the packet at the head of the queue is blocked from moving through the switch fabric because its destination output port is busy, all subsequent packets behind it are prevented from moving forward, even if their requested output ports are idle!

### Q21: What is Longest Prefix Matching? [Slides 19–23]
- When a router looks up a destination IP in its forwarding table, multiple prefix entries may match.
- The router selects the entry with the longest subnet mask (most specific matching bit prefix).
- Implemented in hardware using **Ternary Content Addressable Memory (TCAM)** in a single clock cycle.

### Q22: Why is NAT (Network Address Translation) Controversial? [Slide 89]
1. **Layer Violation**: Routers are Layer 3 devices and should not inspect or modify Layer 4 port numbers.
2. **Violates End-to-End Principle**: Interferes with P2P applications and incoming server connections because internal hosts cannot be directly addressed from outside.
3. **Address Shortage Fix**: Opponents argue NAT delayed the proper adoption of IPv6.

### Q23: Compare IPv4 and IPv6 [Slides 45, 90–93]
| Feature | IPv4 | IPv6 |
|---|---|---|
| **Address Space** | 32 bits ($4.3 \times 10^9$ addresses) | 128 bits ($3.4 \times 10^{38}$ addresses) |
| **Header Length** | Variable (20–60 bytes) | Fixed (40 bytes) |
| **Header Checksum** | Included (recomputed at every hop) | Removed (speeds up router processing) |
| **Fragmentation** | Allowed at routers and source | Only allowed at sending source host |
| **Configuration** | Manual or DHCP | Stateless Address Autoconfiguration (SLAAC) or DHCPv6 |

---

<a name="chapter-5-theory"></a>
## CHAPTER 5: NETWORK LAYER (CONTROL PLANE)

### Q24: Compare Link-State (LS) and Distance Vector (DV) Routing Algorithms [Slides 13, 34–37]
| Dimension | Link-State (Dijkstra) | Distance Vector (Bellman-Ford) |
|---|---|---|
| **Topology Knowledge** | Global knowledge; all nodes have identical link map | Local knowledge; nodes only know costs to neighbors |
| **Message Complexity** | $O(N \cdot E)$ link state broadcast packets | Updates exchanged only between direct neighbors |
| **Convergence Speed** | Fast ($O(N^2)$ algorithm) | Slower; can form routing loops (count-to-infinity) |
| **Robustness** | High; incorrect node calculation is localized | Low; incorrect route advertised poisons whole network |

### Q25: Why are Intra-AS and Inter-AS Routing Protocols Different? [Slide 88]
1. **Policy**:
   - **Inter-AS (BGP)**: Admin policy dominates; an ISP must enforce which traffic can transit its network and hide internal topology.
   - **Intra-AS (OSPF)**: Single administrative domain; policy is not an issue.
2. **Scale**:
   - Hierarchical routing in BGP saves routing table memory across millions of destination prefixes worldwide.
3. **Performance**:
   - Intra-AS focuses purely on shortest/fastest paths.
   - Inter-AS routes are determined by business contracts, peering agreements, and policy rather than shortest path.

### Q26: Why Logically Centralized SDN (Software Defined Networking)? [Slide 93]
1. **Easier Network Management**: Eliminates router-by-router misconfigurations; central controller orchestrates network-wide flows.
2. **Flexible Traffic Engineering**: Can split traffic across multiple paths (load balancing) or route different applications along distinct paths, which is impossible with traditional destination-based routing.
3. **Open Programmability**: Decouples hardware switches (commodity OpenFlow data plane) from software control algorithms.

---

<a name="chapter-6-theory"></a>
## CHAPTER 6: LINK LAYER & LANS

### Q27: What is the difference between a MAC Address and an IP Address? [Slides 43–45]
- **MAC Address (48 bits)**:
  - Flat, hardware address burned into NIC ROM.
  - Portable: remains identical regardless of where the device connects.
  - Used locally to transfer frames between interfaces on the **same subnet**.
- **IP Address (32 bits)**:
  - Hierarchical logical address assigned based on network attachment.
  - Non-portable: changes when host moves to another subnet.
  - Used globally to route datagrams across **different networks**.

### Q28: How does ARP (Address Resolution Protocol) work? [Slides 46–49]
1. Host A wants to send a datagram to Host B on the same LAN but only knows B's IP address.
2. Host A checks its local ARP cache. If not found, Host A broadcasts an **ARP Query packet** encapsulated in an Ethernet broadcast frame (Destination MAC `FF-FF-FF-FF-FF-FF`).
3. All hosts receive the frame; only Host B (matching the queried IP) responds with an **ARP Reply** containing its MAC address (sent unicast to A).
4. Host A caches the `<IP, MAC, TTL>` mapping in its ARP table.

### Q29: Compare Hubs, Switches, and Routers [Slide 75]
| Device | OSI Layer | Forwards Based On | Collision Domains | Broadcast Domains |
|---|---|---|---|---|
| **Hub** | Layer 1 (Physical) | Repeats bits blindly | 1 shared collision domain | 1 broadcast domain |
| **Switch**| Layer 2 (Data Link) | Learns & forwards by MAC address| Separates each port into collision domain | 1 broadcast domain |
| **Router**| Layer 3 (Network) | Routing table by IP address | Separate collision domains | Separates broadcast domains |

### Q30: Why VLANs (Virtual LANs)? [Slides 77–80]
- **Traffic Isolation**: Prevents broadcast storms (ARP/DHCP) from flooding the entire enterprise LAN.
- **Security & Privacy**: Isolates sensitive departments (e.g., Finance vs. CS students).
- **Flexible Administration**: Users moving physical offices do not need physical cable re-patching; switch ports can be dynamically assigned to VLANs in software.

---

<a name="chapter-7-theory"></a>
## CHAPTER 7: WIRELESS & MOBILE NETWORKS

### Q31: What causes Wireless Link Imperfections? [Slides 19–26]
1. **Path Loss & Attenuation**: Signal power decreases with distance ($\sim d^2$ to $d^6$).
2. **Multipath Propagation**: Radio waves reflect off walls and ground, causing delayed signals to interfere at the receiver.
3. **Interference from other sources**: 2.4 GHz ISM band is shared by Wi-Fi, Bluetooth, microwave ovens, baby monitors.
4. **Hidden Terminal Problem**: Terminal A and C cannot hear each other due to distance/obstacles, but both transmit simultaneously to B, causing collisions at B.

### Q32: Difference between CSMA/CD and CSMA/CA [Slides 31–35, Ch 7]
- **CSMA/CD (Wired Ethernet)**: Transceiver can listen while transmitting. If collision occurs, signal is immediately aborted and jam signal sent.
- **CSMA/CA (802.11 Wi-Fi)**: Transceiver cannot detect collisions while sending (transmission power drowns out weak received signal). Hence, it uses **Collision Avoidance** via random backoff and optional RTS/CTS (Request-to-Send / Clear-to-Send) handshakes.

### Q33: Compare LEO and GEO Satellites [Slides 146–147]
| Satellite Type | Altitude | One-Way Propagation Delay | Advantages / Disadvantages |
|---|---|---|---|
| **GEO** (Geostationary) | $36,000\text{ km}$ | $\approx 270\text{ ms}$ | Stationary over equator; huge latency; large dish needed |
| **LEO** (Low Earth Orbit)| $500 - 1500\text{ km}$ | $\approx 20 - 40\text{ ms}$ | Low latency (broadband); satellites move fast; requires constellation (Starlink) |

---

<a name="chapter-8-theory"></a>
## CHAPTER 8: NETWORK SECURITY

### Q34: What are the Four Pillars of Network Security? [Slide 4]
1. **Confidentiality**: Only intended sender and receiver can read message contents (achieved via encryption).
2. **Authentication**: Sender and receiver confirm the true identities of each other.
3. **Message Integrity**: Ensures data has not been altered or tampered with in transit.
4. **Access & Availability**: Services remain operational and accessible to legitimate users (mitigating DoS attacks).

### Q35: Compare Symmetric-Key and Public-Key Cryptography [Slides 11, 16–17, 28]
| Feature | Symmetric Cryptography (AES, DES) | Public-Key Cryptography (RSA) |
|---|---|---|
| **Keys** | Single shared secret key | Key pair (Public key + Private key) |
| **Computational Speed**| 100 to 1,000 times faster | Computationally intensive (modular exponentiation) |
| **Key Distribution** | Difficult; secret key must be securely shared | Easy; public key can be published openly |
| **Primary Usage** | Bulk data encryption | Key exchange (Session keys) & Digital signatures |

### Q36: Why does RSA Work and Why is it Secure? [Slides 24, 27]
- **Why it Works**: Based on Euler's theorem in modular arithmetic:
  $$(m^e)^d \equiv m^{ed} \equiv m \pmod n \quad \text{when } ed \equiv 1 \pmod{(p-1)(q-1)}$$
- **Why it is Secure**: Given public key $(n, e)$, breaking RSA requires finding private key $d$, which requires finding $z = (p-1)(q-1)$. Finding $z$ requires factoring the huge modulus $n$ into its prime factors $p$ and $q$, which is computationally intractable for 2048-bit keys.

### Q37: Evolution of Authentication Protocols (Why did ap1.0 to ap5.0 fail?) [Slides 31–41]
- **ap1.0**: Alice says "I am Alice". *(Fails: Trudy simply declares herself Alice)*.
- **ap2.0**: Alice sends IP packet with her source IP. *(Fails: IP spoofing)*.
- **ap3.0**: Alice sends secret password. *(Fails: Playback/replay attack; Trudy records password)*.
- **ap3.1**: Alice sends encrypted password. *(Fails: Playback attack still works!)*.
- **ap4.0**: Bob sends a **Nonce** $R$ (number used once); Alice returns $R$ encrypted with shared secret key. *(Works with symmetric keys, but requires pre-shared secret)*.
- **ap5.0**: Uses nonce with public keys. *(Fails due to Man-in-the-Middle Attack: Trudy intercepts and substitutes her own public key)*.
- **Final Solution**: **Certification Authorities (CAs)** digitally sign public keys to prevent substitution attacks!

---

# PART 3: ALL CODE SNIPPETS & ALGORITHMIC PSEUDOCODE

<a name="code-udp"></a>
## 1. UDP Socket Programming (Python 3) [Slides 117–118]

### UDP Client (`UDPClient.py`):
```python
from socket import *

# Server connection parameters
serverName = 'localhost'   # Or IP address, e.g. '127.0.0.1'
serverPort = 12000

# Create UDP client socket (AF_INET = IPv4, SOCK_DGRAM = UDP)
clientSocket = socket(AF_INET, SOCK_DGRAM)

# Get message from user
message = input('Input lowercase sentence: ')

# Send message to server (must encode string to bytes)
clientSocket.sendto(message.encode(), (serverName, serverPort))

# Receive response from server (buffer size 2048 bytes)
modifiedMessage, serverAddress = clientSocket.recvfrom(2048)

# Print decoded response and close socket
print('From Server:', modifiedMessage.decode())
clientSocket.close()
```

### UDP Server (`UDPServer.py`):
```python
from socket import *

serverPort = 12000

# Create UDP server socket
serverSocket = socket(AF_INET, SOCK_DGRAM)

# Bind socket to local port number
serverSocket.bind(('', serverPort))

print('The UDP server is ready to receive...')

while True:
    # Read from UDP socket with client address
    message, clientAddress = serverSocket.recvfrom(2048)
    
    # Process message (convert to uppercase)
    modifiedMessage = message.decode().upper()
    
    # Send reply back to client address
    serverSocket.sendto(modifiedMessage.encode(), clientAddress)
```

---

<a name="code-tcp"></a>
## 2. TCP Socket Programming (Python 3) [Slides 121–122]

### TCP Client (`TCPClient.py`):
```python
from socket import *

serverName = 'localhost'
serverPort = 12000

# Create TCP client socket (SOCK_STREAM = TCP)
clientSocket = socket(AF_INET, SOCK_STREAM)

# Initiate 3-way handshake connection to server
clientSocket.connect((serverName, serverPort))

sentence = input('Input lowercase sentence: ')

# Send bytes into TCP byte stream
clientSocket.send(sentence.encode())

# Receive up to 1024 bytes from server
modifiedSentence = clientSocket.recv(1024)

print('From Server:', modifiedSentence.decode())

# Close TCP connection
clientSocket.close()
```

### TCP Server (`TCPServer.py`):
```python
from socket import *

serverPort = 12000

# Create welcoming TCP socket
serverSocket = socket(AF_INET, SOCK_STREAM)

# Bind to port
serverSocket.bind(('', serverPort))

# Listen for incoming connection requests (queue up to 1 connection)
serverSocket.listen(1)

print('The TCP server is ready to receive...')

while True:
    # Accept incoming connection (creates a dedicated connection socket)
    connectionSocket, addr = serverSocket.accept()
    print(f'Connection established with {addr}')
    
    # Receive data from client
    sentence = connectionSocket.recv(1024).decode()
    
    # Capitalize sentence
    capitalizedSentence = sentence.upper()
    
    # Send capitalized string back through connection socket
    connectionSocket.send(capitalizedSentence.encode())
    
    # Close connection socket for this client
    connectionSocket.close()
```

---

<a name="code-socket-events"></a>
## 3. Handling Socket Timeouts [Slide 126]
```python
import socket

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)

# Set timeout of 5.0 seconds on blocking socket operations
s.settimeout(5.0)

try:
    s.connect(('192.0.2.1', 80))
    data = s.recv(1024)
except socket.timeout:
    print('Socket timed out waiting for server response!')
except socket.error as e:
    print(f'Socket connection error: {e}')
finally:
    s.close()
```

---

<a name="code-dijkstra"></a>
## 4. Dijkstra's Link-State Algorithm Pseudocode [Slide 14]
```python
def dijkstra(nodes, source, cost_matrix):
    # 1. Initialization
    N_prime = {source}
    D = {}
    p = {}
    
    for v in nodes:
        if v != source:
            D[v] = cost_matrix.get((source, v), float('inf'))
            p[v] = source if D[v] < float('inf') else None
    
    # 2. Loop
    while len(N_prime) < len(nodes):
        # Find w not in N_prime with minimum D(w)
        w = min((node for node in nodes if node not in N_prime), key=lambda x: D[x])
        N_prime.add(w)
        
        # Update D(v) for all neighbors of w not in N_prime
        for v in nodes:
            if v not in N_prime and (w, v) in cost_matrix:
                if D[w] + cost_matrix[(w, v)] < D[v]:
                    D[v] = D[w] + cost_matrix[(w, v)]
                    p[v] = w
                    
    return D, p
```

---

<a name="code-bellman-ford"></a>
## 5. Bellman-Ford Distance Vector Algorithm Pseudocode [Slide 34]
```python
def distance_vector_update(node_x, neighbors, cost_to_neighbor, neighbor_DVs):
    '''
    Bellman-Ford Equation:
    D_x(y) = min_v { c(x,v) + D_v(y) } for each neighbor v
    '''
    new_Dx = {}
    all_destinations = neighbor_DVs[neighbors[0]].keys()
    
    for y in all_destinations:
        if y == node_x:
            new_Dx[y] = 0
        else:
            new_Dx[y] = min(cost_to_neighbor[v] + neighbor_DVs[v][y] for v in neighbors)
            
    return new_Dx
```

---

<a name="code-csmacd"></a>
## 6. CSMA/CD Binary Exponential Backoff Algorithm [Slide 34]
```python
import random
import time

def csma_cd_backoff(collision_count):
    # m is capped at 10 (IEEE 802.3 Ethernet standard)
    m = min(collision_count, 10)
    
    if collision_count > 16:
        raise Exception("Too many collisions, abort transmission!")
        
    # Choose random slot K from {0, 1, 2, ..., 2^m - 1}
    K = random.randint(0, 2**m - 1)
    
    # Wait time is K * 512 bit times (for 10 Mbps Ethernet, 512 bit times = 51.2 microseconds)
    wait_time_microseconds = K * 51.2
    return wait_time_microseconds
```

---
*End of Comprehensive Study Guide.*
