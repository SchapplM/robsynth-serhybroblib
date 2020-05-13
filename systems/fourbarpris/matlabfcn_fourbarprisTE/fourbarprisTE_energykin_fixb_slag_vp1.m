% Calculate kinetic energy for
% fourbarprisTE
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[GK,GP,HP]';
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
% Datum: 2020-05-07 09:01
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbarprisTE_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisTE_energykin_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisTE_energykin_fixb_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisTE_energykin_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisTE_energykin_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisTE_energykin_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisTE_energykin_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:01:12
% EndTime: 2020-05-07 09:01:13
% DurationCPUTime: 0.11s
% Computational Cost: add. (247->29), mult. (194->55), div. (38->7), fcn. (16->2), ass. (0->22)
t208 = qJ(1) + pkin(3);
t223 = -pkin(2) + t208;
t224 = -pkin(2) - t208;
t234 = (pkin(1) - t223) * (pkin(1) + t224) * (pkin(1) - t224) * (pkin(1) + t223);
t232 = pkin(1) ^ 2;
t216 = sqrt(-t234);
t229 = -t216 / 0.2e1;
t221 = -qJ(1) ^ 2 + (-2 * qJ(1) - pkin(3)) * pkin(3);
t226 = (-pkin(2) ^ 2 + t232);
t195 = t221 + t226;
t217 = (t208 ^ 2);
t206 = 1 / t217;
t228 = t195 * t206;
t225 = qJD(1) / pkin(1);
t205 = 1 / t208;
t194 = -t221 + t226;
t192 = 0.1e1 / t216;
t190 = -rSges(2,1) * t216 - (rSges(2,2) * t194);
t189 = -rSges(2,1) * t194 + rSges(2,2) * t216;
t188 = (-(t194 * t205) / 0.2e1 + (qJ(1) / 0.2e1 - t192 * ((rSges(3,1) * t194) - rSges(3,3) * t216) / 0.2e1) * t228) * t225;
t187 = (t205 * t229 + (rSges(3,1) * t229 + (-qJ(1) / 0.2e1 - rSges(3,3) / 0.2e1) * t194) * t192 * t228) * t225;
t1 = m(3) * (t187 ^ 2 + t188 ^ 2) / 0.2e1 + (-(2 * ((rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3)) * t217) + (m(2) * (-t190 ^ 2 / 0.4e1 - t189 ^ 2 / 0.4e1) / t232 / (t217 ^ 2) / 0.2e1 - ((Icges(3,2) + Icges(2,3)) * t206) / 0.2e1) * (t195 ^ 2)) * qJD(1) ^ 2 / t234;
T = t1;
