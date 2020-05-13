% Calculate kinetic energy for
% fourbarprisDE2
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
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbarprisDE2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_energykin_fixb_slag_vp1: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE2_energykin_fixb_slag_vp1: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_energykin_fixb_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE2_energykin_fixb_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'fourbarprisDE2_energykin_fixb_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'fourbarprisDE2_energykin_fixb_slag_vp1: Icges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:44:56
% EndTime: 2020-05-07 09:44:57
% DurationCPUTime: 0.17s
% Computational Cost: add. (225->53), mult. (238->80), div. (38->8), fcn. (8->2), ass. (0->36)
t250 = qJ(1) + pkin(3);
t276 = -pkin(2) + t250;
t277 = -pkin(2) - t250;
t274 = (0.1e1 / (pkin(1) - t277) / (pkin(1) + t276) / (pkin(1) - t276) / (pkin(1) + t277));
t259 = (pkin(3) ^ 2);
t287 = 3 * t259;
t269 = (t250 ^ 2);
t286 = -2 * t269;
t285 = 2 * qJ(1);
t264 = t259 ^ 2;
t284 = t264 / 0.2e1;
t283 = pkin(3) * qJ(1);
t256 = qJ(1) ^ 2;
t273 = -t256 - t259 - 2 * t283;
t260 = pkin(2) ^ 2;
t261 = pkin(1) ^ 2;
t278 = (t261 - t260);
t236 = t273 + t278;
t282 = rSges(3,1) * t236;
t251 = pkin(1) + pkin(2);
t252 = pkin(1) - pkin(2);
t279 = t252 ^ 2 * t251 ^ 2;
t245 = (-t260 - t261);
t275 = t261 / 0.2e1 - t260 / 0.2e1;
t272 = (1 / t269 ^ 2 / pkin(1) ^ 2 * t274) / 0.2e1;
t267 = t256 ^ 2;
t258 = pkin(3) * t259;
t255 = qJ(1) * t256;
t253 = rSges(3,3) / 0.2e1;
t235 = -t273 + t278;
t233 = sqrt(-0.1e1 / t274);
t232 = -t233 * rSges(2,1) - (rSges(2,2) * t235);
t231 = -rSges(2,1) * t235 + t233 * rSges(2,2);
t230 = ((t235 * t282) / 0.2e1 + (t255 + (0.5e1 / 0.2e1 * pkin(3) + t253) * t256 + (t259 * t285) + t258 / 0.2e1 + t275 * pkin(3) + (t283 + t259 / 0.2e1 - t275) * rSges(3,3)) * t233) * qJD(1);
t229 = (-t233 * t282 / 0.2e1 + t253 * t267 + ((8 * t259 + t245) * t255) + ((t256 * t287) + t284 - t279 / 0.2e1) * rSges(3,3) + (t285 * rSges(3,3) + t245 + (7 * t256)) * t258 + ((t245 * t287 + 3 * t264 + t267) * qJ(1)) + (0.9e1 / 0.2e1 * t267 + (3 * t245 * t256) + t284 + t279 / 0.2e1 + 0.2e1 * rSges(3,3) * t255) * pkin(3)) * qJD(1);
t1 = m(3) * (-t229 ^ 2 - t230 ^ 2) * t272 + ((((rSges(4,1) ^ 2 + rSges(4,2) ^ 2) * m(4) + Icges(4,3)) * t286 * t274) + (((Icges(3,2) + Icges(2,3)) * t274 / t286) + m(2) * (-t232 ^ 2 / 0.4e1 - t231 ^ 2 / 0.4e1) * t272) * (t236 ^ 2)) * qJD(1) ^ 2;
T = t1;
