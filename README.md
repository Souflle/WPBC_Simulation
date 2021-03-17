# WPBC_Simulation
Simulate Wirelessly Powered Backscattering Communication (WPBC) system
* Research for WPBC System for low powered IoT Device

## WPBC SYSTEM

### Configuration
* Reader : Supply Energy & Receive Data from IoT Device
* Tag : Harvest energy from Reader's RF signal and Transfer information to Reader

  ![image](https://user-images.githubusercontent.com/55352279/111420454-9603aa00-872e-11eb-9db0-274a200acc08.png)

### Constraint
* For IoT Device, the power is limited for Sustainable Network.
* Therefore, Information transfer is limited to using passive element ( not using oscilator ) 
* By backscattering communication, Tag can be implemented by passive element. However, in this way, Tag can't estimate channel state

### Requirement
* TX Energy Beamforming for maximize Energy Harvesting rate -> depends on forward channel (F-CSI)
* RX beamforming for reliable communication -> depends on backward channel (B-CSI)


## SYSTEM MODEL

### Channel Model
* Long-term block fading model : Fading coefficient is constant in coherence time
* Path Loss : PL = Ae / (4pi*d^2)

### Slot structure
* CE SLOT (Channel Estimation) : MMSE Estimation of BS-CSI then roghly estimate F-CSI & B-CSI
* WET/IT SLOT : Using Estimated CSI, Tag harvest energy & Information Transfer
 


 





