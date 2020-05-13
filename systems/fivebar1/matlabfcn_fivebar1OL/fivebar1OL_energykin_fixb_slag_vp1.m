% Calculate kinetic energy for
% fivebar1OL
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[AB,AE,BC,CD,ED]';
% m [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fivebar1OL_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_energykin_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_energykin_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'fivebar1OL_energykin_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'fivebar1OL_energykin_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:12:56
% EndTime: 2020-04-27 06:12:56
% DurationCPUTime: 0.04s
% Computational Cost: add. (37->27), mult. (64->54), div. (0->0), fcn. (20->8), ass. (0->23)
t122 = pkin(2) * qJD(1);
t121 = pkin(3) * qJD(3);
t118 = cos(qJ(1));
t117 = cos(qJ(3));
t116 = sin(qJ(1));
t115 = sin(qJ(3));
t114 = qJ(1) + qJ(2);
t113 = qJ(3) + qJ(4);
t112 = qJD(1) + qJD(2);
t111 = qJD(3) + qJD(4);
t110 = cos(t114);
t109 = cos(t113);
t108 = sin(t114);
t107 = sin(t113);
t106 = t118 * rSges(2,1) - t116 * rSges(2,2);
t105 = t117 * rSges(4,1) - t115 * rSges(4,2);
t104 = t116 * rSges(2,1) + t118 * rSges(2,2);
t103 = t115 * rSges(4,1) + t117 * rSges(4,2);
t102 = t118 * t122 + t112 * (-t110 * rSges(3,1) + t108 * rSges(3,2));
t101 = -t116 * t122 - t112 * (-t108 * rSges(3,1) - t110 * rSges(3,2));
t100 = t117 * t121 + t111 * (t109 * rSges(5,1) - t107 * rSges(5,2));
t99 = -t115 * t121 - t111 * (t107 * rSges(5,1) + t109 * rSges(5,2));
t1 = m(3) * (t101 ^ 2 + t102 ^ 2) / 0.2e1 + t112 ^ 2 * Icges(3,3) / 0.2e1 + m(5) * (t100 ^ 2 + t99 ^ 2) / 0.2e1 + t111 ^ 2 * Icges(5,3) / 0.2e1 + (m(2) * (t104 ^ 2 + t106 ^ 2) / 0.2e1 + Icges(2,3) / 0.2e1) * qJD(1) ^ 2 + (m(4) * (t103 ^ 2 + t105 ^ 2) / 0.2e1 + Icges(4,3) / 0.2e1) * qJD(3) ^ 2;
T = t1;
