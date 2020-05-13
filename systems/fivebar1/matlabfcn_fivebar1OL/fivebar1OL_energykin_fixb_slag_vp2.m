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
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-04-27 06:13
% Revision: 970157063c8b7bcb25458a457551a083b21abdbd (2020-04-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = fivebar1OL_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'fivebar1OL_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'fivebar1OL_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'fivebar1OL_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'fivebar1OL_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'fivebar1OL_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'fivebar1OL_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-04-27 06:12:56
% EndTime: 2020-04-27 06:12:56
% DurationCPUTime: 0.04s
% Computational Cost: add. (17->13), mult. (52->32), div. (0->0), fcn. (8->4), ass. (0->7)
t57 = cos(qJ(2));
t56 = cos(qJ(4));
t55 = sin(qJ(2));
t54 = sin(qJ(4));
t53 = qJD(1) + qJD(2);
t52 = qJD(3) + qJD(4);
t1 = (Ifges(2,3) / 0.2e1 + m(3) * (t55 ^ 2 + t57 ^ 2) * pkin(2) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (Ifges(4,3) / 0.2e1 + m(5) * (t54 ^ 2 + t56 ^ 2) * pkin(3) ^ 2 / 0.2e1) * qJD(3) ^ 2 + (Ifges(3,3) * t53 / 0.2e1 + (-mrSges(3,1) * t57 + mrSges(3,2) * t55) * qJD(1) * pkin(2)) * t53 + (Ifges(5,3) * t52 / 0.2e1 + (mrSges(5,1) * t56 - mrSges(5,2) * t54) * qJD(3) * pkin(3)) * t52;
T = t1;
