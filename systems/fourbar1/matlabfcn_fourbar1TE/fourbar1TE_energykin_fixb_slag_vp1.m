% Calculate kinetic energy for
% fourbar1TE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[l1,l2,l3,l4]';
% m [4x1]
%   mass of all robot links (including the base)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-24 19:49
% Revision: 381275ccc1b3107bade6cbef2523183900cd8455 (2020-04-24)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbar1TE_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbar1TE_energykin_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbar1TE_energykin_fixb_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'fourbar1TE_energykin_fixb_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbar1TE_energykin_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbar1TE_energykin_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbar1TE_energykin_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-24 19:48:57
% EndTime: 2020-04-24 19:48:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (308->39), mult. (441->74), div. (18->7), fcn. (110->4), ass. (0->34)
t222 = pkin(4) ^ 2;
t199 = cos(qJ(1));
t216 = pkin(2) * t199;
t214 = (-0.2e1 * t216 + pkin(1)) * pkin(1);
t219 = -pkin(3) + pkin(4);
t220 = -pkin(3) - pkin(4);
t221 = ((pkin(2) - t220) * (pkin(2) + t220) + t214) * ((pkin(2) - t219) * (pkin(2) + t219) + t214);
t198 = sin(qJ(1));
t217 = pkin(2) * t198;
t215 = pkin(2) * qJD(1);
t213 = pkin(3) ^ 2 - t222;
t205 = pkin(2) ^ 2;
t196 = t205 + t214;
t191 = t196 - t213;
t192 = 0.1e1 / t196;
t207 = sqrt(-t221);
t211 = -pkin(1) + t216;
t176 = (0.1e1 - (t191 * t217 - t211 * t207) / t207 * pkin(1) * t192) * qJD(1);
t212 = t176 * t192 / pkin(3);
t208 = t196 ^ 2;
t200 = qJD(1) ^ 2;
t195 = t199 * rSges(2,1) - t198 * rSges(2,2);
t194 = t198 * rSges(2,1) + t199 * rSges(2,2);
t190 = t196 + t213;
t189 = -t211 * rSges(3,1) + rSges(3,2) * t217;
t188 = rSges(3,1) * t217 + t211 * rSges(3,2);
t187 = rSges(4,1) * t217 + t211 * rSges(4,2);
t186 = -t211 * rSges(4,1) + rSges(4,2) * t217;
t180 = (pkin(1) * t199 - pkin(2)) * t207 - pkin(1) * t198 * t190;
t178 = -t191 * t186 + t187 * t207;
t177 = t186 * t207 + t191 * t187;
t175 = t199 * t215 + (t188 * t207 + t190 * t189) * t212 / 0.2e1;
t174 = -t198 * t215 - (-t190 * t188 + t189 * t207) * t212 / 0.2e1;
t1 = m(3) * (t174 ^ 2 + t175 ^ 2) / 0.2e1 + t176 ^ 2 * Icges(3,3) / 0.2e1 + (m(2) * (t194 ^ 2 + t195 ^ 2) + Icges(2,3)) * t200 / 0.2e1 + (m(4) * (-t177 ^ 2 / 0.4e1 - t178 ^ 2 / 0.4e1) / t222 / t208 ^ 2 / 0.2e1 - 0.1e1 / t208 * Icges(4,3) / 0.2e1) * t180 ^ 2 * t200 * t205 / t221;
T = t1;
