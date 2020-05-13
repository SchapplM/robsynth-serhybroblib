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
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-05-07 09:45
% Revision: 70b73dc947b5bf54b9f851309d04479b7d95fc8d (2020-05-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fourbarprisDE2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'fourbarprisDE2_energykin_fixb_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'fourbarprisDE2_energykin_fixb_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'fourbarprisDE2_energykin_fixb_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'fourbarprisDE2_energykin_fixb_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'fourbarprisDE2_energykin_fixb_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'fourbarprisDE2_energykin_fixb_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-05-07 09:44:51
% EndTime: 2020-05-07 09:44:51
% DurationCPUTime: 0.06s
% Computational Cost: add. (105->18), mult. (70->21), div. (25->7), fcn. (2->2), ass. (0->8)
t186 = (qJ(1) + pkin(3));
t190 = -pkin(2) + t186;
t191 = -pkin(2) - t186;
t197 = (0.1e1 / (pkin(1) + t191) / (pkin(1) - t190) / (pkin(1) - t191) / (pkin(1) + t190));
t189 = (t186 ^ 2);
t188 = (qJ(1) ^ 2);
t173 = pkin(1) ^ 2 - pkin(2) ^ 2 - t188 - (2 * qJ(1) + pkin(3)) * pkin(3);
t1 = (m(3) / 0.2e1 + (t173 / t186 * (-1 / t197) ^ (-0.1e1 / 0.2e1) * mrSges(3,1)) + (-(2 * Ifges(4,3) * t189) + (-Ifges(2,3) / 0.2e1 - Ifges(3,2) / 0.2e1 - m(3) * t188 / 0.2e1 - (qJ(1) * mrSges(3,3))) * (t173 ^ 2) / t189) * t197) * qJD(1) ^ 2;
T = t1;
